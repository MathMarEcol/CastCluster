#' CAST algorithm for a given number of clusters.
#'
#' Rather than fixing the threshold and
#' finding k that meets the threshold,
#' this function fixes k regardless of
#' threshold.
#'
#'
#' I tried just allocating all the seeds at once
#' but it tends to find poor local optima
#'
#' @param sim_mat symmetric matrix with values in the range [0,1]
#' @param k number of clusters to fit
#' @param max_iter don't shuffle objects forever
cast_k <- function(sim_mat, k,
                   max_iter = nrow(sim_mat)*2){

  assertthat::assert_that(isSymmetric(unname(sim_mat)))

  assertthat::assert_that(k > 1)
  assertthat::assert_that(k <= nrow(sim_mat))

  aff_thres <- mean(sim_mat) ##need a stopping point. Subdivide if too broad


  ##get an initial CAST clustering and stabilize
  cast_ob <- cast_alg(sim_mat, aff_thres)
  cast_ob <- cast_stabilize(cast_ob, aff_thres = 0, sim_mat)

  ##Add or remove clusters until k is reached.
  while(length(cast_ob) != k){
    clust_sim <- aff_cluster_between(cast_ob, sim_mat)
    clust_sim_mat <- matrix(clust_sim$affs, nrow = length(cast_ob), ncol =length(cast_ob) )
    if(length(cast_ob) > k){
      ## ##Merge two most similar clusters
      ## clust_sim <- clust_sim[!(clust_sim$Var1 == clust_sim$Var2),]
      ## clust_merge <- clust_sim[which.max(clust_sim$affs),]
      ## cast_ob[[clust_merge$Var1]] <- c(cast_ob[[clust_merge$Var1]], cast_ob[[clust_merge$Var2]])
      ## cast_ob[[clust_merge$Var2]] <- NULL
      ## Remove weakest cluster
      clust_sim <- clust_sim[(clust_sim$Var1 == clust_sim$Var2),]
      clust_remove <- which.min(clust_sim$affs)

      ##I need to follow the logic of cast_compact
      ##Simplify, I just pick one move at a time until all sites are gone
      ##then use cast_stabilize() to finish then job

      ## matrix, sum affinity and cluster N matrix
      clust_aff_sum <- aff_clust_sum(sim_mat, cast_ob)
      clust_aff_n <- do.call(c, lapply(cast_ob, length))
      clust_aff <- clust_aff_sum * rep(1 / clust_aff_n, each = nrow(clust_aff_sum))
      clust_aff[, clust_remove] <- -2

      ##get mapping from elements to clusters
      clust_ind <- lapply(seq_along(cast_test), function(clust, cast_test){
        if(length(cast_test[[clust]]) > 0) {
          return(data.frame(elem = cast_test[[clust]], clust = clust))
        } else {
          return(NULL)
        }
      }, cast_test = cast_test)
      clust_ind <- do.call(rbind, clust_ind)
      clust_ind <- clust_ind[order(clust_ind$elem), ]
      while(length(cast_ob[[clust_remove]]) != 0 ){
        ind <- clust_ind$elem + (clust_ind$clust - 1) * nrow(clust_ind)
        clust_current <- clust_aff[ind]
        gain <- clust_aff - clust_current
        max_ind <- which.max(gain)
        max_ind_zero <- (max_ind - 1)
        i <- 1 + (max_ind_zero %% nrow(clust_ind)) #modulo is zero based indexing, R is 1 based indexing
        from <- clust_ind$clust[i]
        assertthat::assert_that(from == clust_remove)
        to <- (max_ind - i)/nrow(clust_ind) + 1
        upd <- data.frame(i = i, from = from, to = to, gain = gain[max_ind])


        cast_ob[[upd$to]] <- c(cast_ob[[upd$to]], upd$i)
        cast_ob[[upd$from]] <- cast_ob[[ upd$from ]][cast_ob[[ upd$from ]] != upd$i  ] 

        clust_ind[upd$i, 2] <- upd$to

        clust_aff_n[upd$to] <- clust_aff_n[upd$to] + 1
        clust_aff_sum[, upd$to] <- clust_aff_sum[, upd$to] + sim_mat[, upd$i]
        clust_aff[, upd$to] <- clust_aff_sum[, upd$to] * (1 / clust_aff_n[upd$to])

      }


      ##Stabilise
      cast_ob <- cast_stabilize(cast_ob, aff_thres = 0, sim_mat)
    } else {
      ##Seed the weakest object into a new cluster
      affs <- castcluster::aff_clust_all(sim_mat, cast_ob)
      new_seed <- which.min(affs)
      new_seed_zero <- (new_seed - 1)
      i <- 1 + (new_seed_zero %% nrow(affs)) #modulo is zero based indexing, R is 1 based indexing
      from <- (max_ind - i)/nrow(affs) + 1
      cast_ob[[from]] <- cast_ob[[from]][cast_ob[[from]] != i]
      cast_ob[[length(cast_ob) + 1]] <- c(i)

      ##Stabilise
      cast_ob <- cast_stabilize(cast_ob, aff_thres = 0, sim_mat)

    }

  }

  return(cast_ob)

}
