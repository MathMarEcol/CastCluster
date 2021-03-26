#' Stabilize cluster membership
#'
#' cast_stabilize() works iteratively, updating one site at a time
#' and recalculating the affinities each time.
#' cast_stabilize() will only make an update if it either maintains
#' all the average within cluster affinities above the threshold, or
#' increases the within cluster affinities if one is already below the
#' threshold.
#'
#' cast_stabilize_batch() does a bulk update, finding all sites that
#' could be assigned to another cluster before recalulating the affinities.
#' cast_stabilize_batch() always proceeds without checking the impact
#' the update has on within-cluster affinities.
cast_stabilize_batch <- function(cast_obj, aff_thres, sim_mat, max_iter = 20){
  iter <- 1

  while(iter <= max_iter){
    ##For each vertex, find affinity to other clusters
    ##Cast_obj is not huge, but constantly recreating it when I am just
    ##flipping bools seems wasteful
    ##This approach tests each site, and updates clustering at end
    updates <- lapply(seq_along(cast_obj[[1]]), function(u, cast_obj, sim_mat){
      clust_id <- which(sapply(cast_obj, function(clust, u){
        u %in% clust
      }, u = u))
      assertthat::assert_that(length(clust_id) == 1)
      ##u belongs to clust_id
      clust_aff <- lapply(cast_obj, function(clust, u, sim_mat){
          ##get the affinity to each cluster
          if(any(clust)){
            return(mean(sim_mat[u, clust]))
          } else {
            ##empty clusters have 0 affinity
            return(0)
          }
        }, u = u, sim_mat = sim_mat)
        if(which.max(clust_aff) != clust_id){
          return(data.frame(u = u, old = clust_id, new = which.max(clust_aff)))
        } else {
          return(NULL)
        }
      }, cast_obj = cast_obj, sim_mat = sim_mat)
    updates <- do.call("rbind", updates)
    ##Apply updates
    if(!is.null(updates)){
      message("iteration [", iter, "] reassigned [", nrow(updates),"] samples")
      for(upd in 1:nrow(updates) ){
        cast_obj[[ updates[upd,"old"] ]] <- cast_obj[[ updates[upd,"old"] ]][cast_obj[[ updates[upd,"old"] ]] != updates[upd,"u"]  ] 
        cast_obj[[updates[upd,"new"]]] <- c(cast_obj[[updates[upd,"new"]]], updates[upd,"u"])
      }
    } else {
      break
    }
    iter <- iter + 1
  }
  return(cast_obj)
}
cast_stabilize <- function(cast_obj, aff_thres, sim_mat, max_iter = nrow(sim_mat)*2){

  sim_mat_diag <- diag(sim_mat)
  diag(sim_mat) <- NA
  ##matrix, mean affinity of each site to each cluster
  ##two cols must be updated each iteration
  ## ## clust_aff <- aff_clust_all(sim_mat, cast_test)
  ## clust_aff_sum <- aff_clust_sum(sim_mat, cast_obj)
  ## ## matrix, sum affinity and cluster N matrix
  ## clust_aff_n <- do.call(c, lapply(cast_obj, length))
  ## clust_aff <- clust_aff_sum * rep(1 / clust_aff_n, each = nrow(clust_aff_sum))
  clust_aff <- castcluster:::.aff_clust_mean(sim_mat, cast_obj, aff_thres)

  ##maps element to cluster
  ##can be updated with a single value change each loop
  clust_ind <- lapply(seq_along(cast_obj), function(clust, cast_obj){
    data.frame(elem = cast_obj[[clust]], clust = clust)
  }, cast_obj = cast_obj)
  clust_ind <- do.call(rbind, clust_ind)
  clust_ind <- clust_ind[order(clust_ind$elem), ]

  ##sites that should not be moved, to preserve aff_thres
  locked_sites <- c()
  within_clust_affs <- do.call(c, lapply(seq_along(cast_obj), function(clust, cast_obj, clust_aff){
    mean(clust_aff[ cast_obj[[clust]], clust])
  }, clust_aff = clust_aff, cast_obj = cast_obj))

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
    max_ind_zero <- (max_ind -1)
    i <- 1 + (max_ind_zero %% nrow(clust_ind)) #modulo is zero based indexing, R is 1 based indexing
    from <- clust_ind$clust[i]
    to <- (max_ind - i)/nrow(clust_ind) + 1
    upd <- data.frame(i = i, from = from, to = to, gain = gain[max_ind])

    ##Apply updates
    if(upd$gain > 0){
      ## message("applying update [", iter,"]")
      ## message(paste(names(upd), collapse =" - "))
      ## message(paste(upd, collapse = " - "))


      to_clust <- c(cast_obj[[upd$to]], upd$i)
      from_clust <- cast_obj[[ upd$from ]][cast_obj[[ upd$from ]] != upd$i  ] 
      ## clust_sum_to <- clust_aff_sum[, upd$to] + sim_mat[, upd$i]
      ## clust_aff_to <- 1/(clust_aff_n[upd$to] + 1) *  mean(clust_sum_to[to_clust])
      clust_aff_to <- mean(castcluster:::.aff_clust_mean(sim_mat, list(to_clust), aff_thres)[to_clust])
      ## clust_sum_from <- clust_aff_sum[, upd$from] - sim_mat[, upd$i]
      ## clust_aff_from <- 1/(clust_aff_n[upd$from] - 1) * mean(clust_sum_from[from_clust])
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

      is_update <- (is_above_to || is_improve_to)  && (is_above_from || is_improve_from)


      ##Update is also good if the affinity of clusters started below threshold
      ##and the update has not made the situation worse.

      if(is_update){
        cast_obj[[upd$to]] <- to_clust
        cast_obj[[upd$from]] <- from_clust

        clust_ind[upd$i, 2] <- upd$to

        ## clust_aff_n[upd$to] <- clust_aff_n[upd$to] + 1
        ## clust_aff_sum[, upd$to] <- clust_sum_to
        ## clust_aff[, upd$to] <- clust_aff_sum[, upd$to] * (1 / clust_aff_n[upd$to])
        clust_aff[, upd$to] <- castcluster:::.aff_clust_mean(sim_mat, cast_obj[upd$to], aff_thres)

        ## clust_aff_n[upd$from] <- clust_aff_n[upd$from] - 1
        ## clust_aff_sum[, upd$from] <- clust_sum_from
        ## clust_aff[, upd$from] <- clust_aff_sum[, upd$from] * (1 / clust_aff_n[upd$from])
        clust_aff[, upd$from] <- castcluster:::.aff_clust_mean(sim_mat, cast_obj[upd$from], aff_thres)

        within_clust_affs[upd$from] <- clust_aff_from
        within_clust_affs[upd$to] <- clust_aff_to

        locked_sites <- c()

        if(length(from_clust) == 0) {
          ##from cluster is empty: remove
          message("cluster [", upd$from, "] has been removed.")
          cast_obj[[upd$from]] <- NULL
          clust_aff <- clust_aff[, -upd$from]
          within_clust_affs  <- within_clust_affs[-upd$from]

          clust_ind <- lapply(seq_along(cast_obj), function(clust, cast_obj){
            data.frame(elem = cast_obj[[clust]], clust = clust)
          }, cast_obj = cast_obj)
          clust_ind <- do.call(rbind, clust_ind)
          clust_ind <- clust_ind[order(clust_ind$elem), ]

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

    } else {
      break
    }
    iter <- iter + 1
  }
  return(cast_obj)
}
