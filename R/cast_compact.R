# Copyright 2017-2024 Philip Dyer
# SPDX-License-Identifier: GPL-3.0-only
#' cast_compact takes a stabilized cast object, and deletes
#' one cluster, assignes the sites to nearest clusters,
#' stabilizes again, and
#' checks whether the average cluster affinity
#' remains above the threshold. If the threshold is still valid,
#' keep the changed cast algorithm, otherwise revert and try another cluster.
cast_compact_batch <- function(cast_ob, sim_mat, aff_thres, max_iter = nrow(sim_mat)*2){
  ##sort by size
  cluster_size <- data.frame(clust = 1:length(cast_ob), size = sapply(cast_ob, length))
  cluster_size <- cluster_size[order(cluster_size$size),]
  i <- 1
  repeat{
    if(length(cast_ob) <= 1){
      ##Stop at one cluster
      break
    }
    cast_test <- cast_ob
    ##assign affinities if cluster i to other clusters
    clust_id <- cluster_size$clust[i]
    updates <- lapply(seq_along(cast_ob[[clust_id]]), function(u, cast_ob, sim_mat){
      clust_aff <- sapply(cast_ob, function(clust, site, sim_mat){
          ##get the affinity to each cluster
          if(length(clust) > 0){
            return(mean(sim_mat[site, clust]))
          } else {
            ##empty clusters have 0 affinity
            return(0)
          }
      }, site = cast_ob[[clust_id]][u], sim_mat = sim_mat)
      new_clust <- which(clust_aff == max(clust_aff[-clust_id]))[1]
      return(data.frame(site = cast_ob[[clust_id]][u], old = clust_id, new = new_clust))
      }, cast_ob = cast_ob, sim_mat = sim_mat)
    updates <- do.call("rbind", updates)
    ##Apply updates
      #message("reassigned [", nrow(updates),"] samples from cluster [", clust_id, "]")
      for(upd in 1:nrow(updates) ){
        cast_test[[ updates[upd,"old"] ]] <- cast_test[[ updates[upd,"old"] ]][cast_test[[ updates[upd,"old"] ]] != updates[upd,"site"]  ] 
        cast_test[[updates[upd,"new"]]] <- c(cast_test[[updates[upd,"new"]]], updates[upd,"site"])
      }

    cast_test[[clust_id]] <- NULL

    ##test new cluster affinities
    new_aff <- do.call("c", aff_clust_inner(cast_obj = cast_test, sim_mat = sim_mat))
    #message("min_new_aff: [", min(new_aff), "] and threshold of [", aff_thres, "]. attempt: [", i, "]")
    cast_test_stab <- cast_stabilize_batch(cast_obj = cast_test, aff_thres = aff_thres, sim_mat = sim_mat, max_iter = max_iter )
    new_aff_stab <- do.call("c", aff_clust_inner(cast_obj = cast_test_stab, sim_mat = sim_mat))
    #message("min_new_aff_stab: [", min(new_aff_stab), "] and threshold of [", aff_thres, "]")
    if(min(new_aff_stab) >= aff_thres && min(new_aff) >= aff_thres){
      cast_ob <- cast_test_stab
      cluster_size <- data.frame(clust = 1:length(cast_ob), size = sapply(cast_ob, length))
      cluster_size <- cluster_size[order(cluster_size$size),]
      i <- 0
    }

    i <- i + 1
    if(i > nrow(cluster_size)){
      break
    }
  }
  return(cast_ob)
}

#' cast_compact
#'
#' takes a stabilized cast object, and deletes
#' one cluster, assignes the sites to nearest clusters,
#' stabilizes again, and
#' checks whether the average cluster affinity
#' remains above the threshold. If the threshold is still valid,
#' keep the changed cast algorithm, otherwise revert and try another cluster.
#'
#' @export
cast_compact <- function(cast_ob, sim_mat, aff_thres, max_iter = nrow(sim_mat)*2){
  ##sort by size
  clust_remove <- 1
  sim_mat_diag <- diag(sim_mat)
  diag(sim_mat) <- NA


  ##matrix, mean affinity of each site to each cluster
  ##two cols must be updated each iteration
  ## clust_aff <- aff_clust_all(sim_mat, cast_test)
  ## clust_aff_sum_master <- aff_clust_sum(sim_mat, cast_ob)
  ## ## matrix, sum affinity and cluster N matrix
  ## clust_aff_n_master <- do.call(c, lapply(cast_ob, length))
  ## clust_aff_master <- clust_aff_sum_master * rep(1 / clust_aff_n_master, each = nrow(clust_aff_sum_master))
  clust_aff_master <- castcluster:::.aff_clust_mean(sim_mat, cast_ob, aff_thres)
  ## print(range(clust_aff_master))
  repeat{
    if(length(cast_ob) <= 1){
      ##Stop at one cluster
      break
    }
    ##assign affinities if cluster i to other clusters

    cast_test <- cast_ob
    cluster_size <- data.frame(clust = 1:length(cast_test), size = sapply(cast_test, length))
    cluster_size <- cluster_size[order(cluster_size$size),]
    clust_id <- cluster_size$clust[clust_remove]
    is_removed <- FALSE



    ##maps element to cluster
    ##can be updated with a single value change each loop
    clust_ind <- lapply(seq_along(cast_test), function(clust, cast_test){
      if(length(cast_test[[clust]]) > 0) {
        return(data.frame(elem = cast_test[[clust]], clust = clust))
      } else {
        return(NULL)
      }
    }, cast_test = cast_test)
    clust_ind <- do.call(rbind, clust_ind)
    clust_ind <- clust_ind[order(clust_ind$elem), ]
    ##best move is difference between current group and other group
    ## Gain can range from -1 to 1, so
    ##using an absolute of -2 means points in clust_id will move before anything else
    ##with a range of 1 to 3
    clust_aff <- clust_aff_master
    clust_aff[, clust_id] <- -2

    ## print(range(clust_aff))
    ##sites that should not be moved, to preserve aff_thres
    locked_sites <- numeric(0)
    within_clust_affs <- do.call(c, lapply(seq_along(cast_test), function(clust, cast_ob, clust_aff){
    mean(clust_aff[cast_ob[[clust]], clust])
  }, clust_aff = clust_aff, cast_ob = cast_test))

    ##disband and stabilize
    iter <- 1
    while(iter <= max_iter){
      ##single step updates require some conserved data.

      ##find maximal out of cluster affinity

      ##only have to update clusters that change

      ##candidates are sites from other clusters that exceed threshold

      ##data.frame, element with highest affinity gain by changing clusters
      ##must be recalculated every loop
      ind <- clust_ind$elem + (clust_ind$clust - 1) * nrow(clust_ind)
      clust_current <- clust_aff[ind]
      gain <- clust_aff - clust_current
      if(length(locked_sites ) > 0) {
        gain[locked_sites, ] <- 0
      }
      max_ind <- which.max(gain)
      max_ind_zero <- (max_ind - 1)
      i <- 1 + (max_ind_zero %% nrow(clust_ind)) #modulo is zero based indexing, R is 1 based indexing
      from <- clust_ind$clust[i]
      to <- (max_ind - i)/nrow(clust_ind) + 1
      upd <- data.frame(i = i, from = from, to = to, gain = gain[max_ind])
      ## print(upd)
      ## print(range(clust_aff))

      ##Apply updates
      if(upd$gain > 0){
        ## message("applying update [", iter,"]")
        ## message(paste(names(upd), collapse =" - "))
        ## message(paste(upd, collapse = " - "))

        to_clust <- c(cast_test[[upd$to]], upd$i)
        from_clust <- cast_test[[ upd$from ]][cast_test[[ upd$from ]] != upd$i  ]
        clust_aff_to <- mean(castcluster:::.aff_clust_mean(sim_mat, list(to_clust), aff_thres)[to_clust])
        clust_aff_from <- mean(castcluster:::.aff_clust_mean(sim_mat, list(from_clust), aff_thres)[from_clust])

        ## Cluster may be eliminated
        if(length(from_clust) == 0){
          clust_aff_from <- aff_thres
        }

        ##Update is good if, after the update, both clusters have either affinity above
        ## aff_thres, or improve affinity.
        is_above_to <- clust_aff_to >= aff_thres
        is_above_from <- clust_aff_from >= aff_thres
        is_improve_to <- within_clust_affs[upd$to] <= clust_aff_to
        is_improve_from <- within_clust_affs[upd$from] <= clust_aff_from

        if(!is_removed && upd$from == clust_id){
          ##Update is moving a site from the disbanded cluster


          ##Update is good if, after the update, the to cluster has affinity above
          is_update <- (is_above_to || is_improve_to)
        } else {
          ##Update is moving a site between two surviving clusters, to stablize

          ##Update is also good if the affinity of clusters started below threshold
          ##and the update has not made the situation worse.
          is_update <- (is_above_to || is_improve_to)  && (is_above_from || is_improve_from)
        }

        if(is_update){
          cast_test[[upd$to]] <- to_clust
          cast_test[[upd$from]] <- from_clust

          clust_ind[upd$i, 2] <- upd$to

          ## clust_aff_n[upd$to] <- clust_aff_n[upd$to] + 1
          ## clust_aff_sum[, upd$to] <- clust_sum_to
          ## clust_aff[, upd$to] <- clust_aff_sum[, upd$to] * (1 / clust_aff_n[upd$to])
          clust_aff[, upd$to] <- castcluster:::.aff_clust_mean(sim_mat, cast_test[upd$to], aff_thres)
          within_clust_affs[upd$to] <- clust_aff_to

          if(is_removed || upd$from != clust_id){
            clust_aff[, upd$from] <- castcluster:::.aff_clust_mean(sim_mat, cast_test[upd$from], aff_thres)
            within_clust_affs[upd$from] <- clust_aff_from
          }
          locked_sites <- numeric(0)

          if(length(from_clust) == 0) {
            ##from cluster is empty: remove
            #message("cluster [", upd$from, "] has been removed.")
            cast_test[[upd$from]] <- NULL
            clust_aff <- clust_aff[, -upd$from]
            within_clust_affs  <- within_clust_affs[-upd$from]

            clust_ind <- lapply(seq_along(cast_test), function(clust, cast_ob){
              data.frame(elem = cast_ob[[clust]], clust = clust)
            }, cast_ob = cast_test)
            clust_ind <- do.call(rbind, clust_ind)
            clust_ind <- clust_ind[order(clust_ind$elem), ]
            if (!is_removed && upd$from == clust_id) {
              ##Target cluster has been removed
              is_removed <- TRUE
            }
          }
        } else {
          if(length(locked_sites ) > 0) {
            locked_sites <- c(locked_sites, upd$i)
          } else {
            locked_sites <- upd$i
          }
          ## message("clust_aff_to [", clust_aff_to ,"] may have dropped below threshold [", aff_thres, "]")
          ## message("clust_aff_from [", clust_aff_from ,"] may have dropped below threshold [", aff_thres, "]")
        }
        ## print(range(clust_aff))
      } else {
        break
      }
      iter <- iter + 1
    }
    ##check that the site has been disbanded
    if(is_removed) {
      ##cluster was successfully disbanded within max_iter
      ## cast_test[clust_id] <- NULL
      cast_ob <- cast_test

      clust_aff_master <- clust_aff
      clust_remove <- 0

      #message("Disbanded cluster [", clust_id, "]")
    } else {
      ##message("Could not disband cluster [", clust_id, "]")
      #message("clust_remove: [", clust_remove, "]")
    }

    ##try next cluster
    clust_remove <- clust_remove + 1
    if(clust_remove > nrow(cluster_size)){
      break
    }
  }
  #message("New number of clusters: [", length(cast_ob), "]")
  return(cast_ob)
}

