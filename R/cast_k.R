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

  n <- nrow(sim_mat)
  spares <- seq.int(n)

  seeds <- sample.int(n, k, replace = FALSE)
  cast_ob <- lapply(seeds, function(x){c(x)})
  spares <- spares[-seeds]

  ##matrix, mean affinity of each site to each cluster
  ##two cols must be updated each iteration
  ## clust_aff <- aff_clust_all(sim_mat, cast_test)
  clust_aff_sum_master <- aff_clust_sum(sim_mat, cast_ob)
  ## matrix, sum affinity and cluster N matrix
  clust_aff_n_master <- do.call(c, lapply(cast_ob, length))
  clust_aff_master <- clust_aff_sum_master * rep(1 / clust_aff_n_master, each = nrow(clust_aff_sum_master))
  clust_aff_sum <- clust_aff_sum_master
  clust_aff_n <- clust_aff_n_master
  clust_aff <- clust_aff_master
  ##Assign all spares, iteratively, greedy
  while(length(spares) > 0) {

    clust_aff_spares <- cbind(seq.int(n), clust_aff)[spares, , drop = FALSE]

    max_ind <- which.max(clust_aff_spares[,-1])
    max_ind_zero <- (max_ind - 1)
    i <- 1 + (max_ind_zero %% nrow(clust_aff_spares)) #modulo is zero based indexing, R is 1 based indexing
    to <- (max_ind - i)/nrow(clust_aff_spares) + 1
    upd <- data.frame(spare = clust_aff_spares[i,1], to = to)

    to_clust <- c(cast_ob[[upd$to]], upd$spare)
    spares <- spares[-i]

    clust_sum_to <- clust_aff_sum[, upd$to] + sim_mat[, upd$spare]
    clust_aff_to <- (1/(clust_aff_n[upd$to] + 1)) *  mean(clust_sum_to[to_clust])

    cast_ob[[upd$to]] <- to_clust

    clust_aff_n[upd$to] <- clust_aff_n[upd$to] + 1
    clust_aff_sum[, upd$to] <- clust_sum_to
    clust_aff[, upd$to] <- clust_aff_sum[, upd$to] * (1 / clust_aff_n[upd$to])
  }



  ##Stabilize the clustering
  clust_ind <- lapply(seq_along(cast_ob), function(clust, cast_ob){
    data.frame(elem = cast_ob[[clust]], clust = clust)
  }, cast_ob = cast_ob)
  clust_ind <- do.call(rbind, clust_ind)
  clust_ind <- clust_ind[order(clust_ind$elem), ]

  iter <- 1
  while(iter <= max_iter){
    ind <- clust_ind$elem + (clust_ind$clust - 1) * nrow(clust_ind)
    clust_current <- clust_aff[ind]
    gain <- clust_aff - clust_current

    max_ind <- which.max(gain)
    max_ind_zero <- (max_ind -1)
    i <- 1 + (max_ind_zero %% nrow(clust_ind)) #modulo is zero based indexing, R is 1 based indexing
    from <- clust_ind$clust[i]
    to <- (max_ind - i)/nrow(clust_ind) + 1
    upd <- data.frame(i = i, from = from, to = to, gain = gain[max_ind])

    if(upd$gain <= 0){
      ##All moves would take an object out of the best cluster.
      break
    } else {

      to_clust <- c(cast_ob[[upd$to]], upd$i)
      from_clust <- cast_ob[[ upd$from ]][cast_ob[[ upd$from ]] != upd$i  ] 
      clust_sum_to <- clust_aff_sum[, upd$to] + sim_mat[, upd$i]
      clust_aff_to <- 1/(clust_aff_n[upd$to] + 1) *  mean(clust_sum_to[to_clust])
      clust_sum_from <- clust_aff_sum[, upd$from] - sim_mat[, upd$i]
      clust_aff_from <- 1/(clust_aff_n[upd$from] - 1) * mean(clust_sum_from[from_clust])

      cast_ob[[upd$to]] <- to_clust
      cast_ob[[upd$from]] <- from_clust

      clust_ind[upd$i, 2] <- upd$to

      clust_aff_n[upd$to] <- clust_aff_n[upd$to] + 1
      clust_aff_sum[, upd$to] <- clust_sum_to
      clust_aff[, upd$to] <- clust_aff_sum[, upd$to] * (1 / clust_aff_n[upd$to])

      clust_aff_n[upd$from] <- clust_aff_n[upd$from] - 1
      clust_aff_sum[, upd$from] <- clust_sum_from
      clust_aff[, upd$from] <- clust_aff_sum[, upd$from] * (1 / clust_aff_n[upd$from])

      all_clust_affs[upd$from] <- clust_aff_from
      all_clust_affs[upd$to] <- clust_aff_to

    }

  }



}
