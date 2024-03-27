# Copyright 2017-2024 Philip Dyer
# SPDX-License-Identifier: GPL-3.0-only
#'
#' @export
cast_alg <- function(sim_mat, aff_thres){
##rather than wrestle with integer sets, I will use boolean vectors

  ##Set the diagonal of sim_mat to 0
  ##combined with division by n-1 during the removal stage
  ##this has the effect of ignoring self affinity
  ## when calculating cluster affinity.
  sim_mat_diag <- diag(sim_mat)
  diag(sim_mat) <- 0

  ##initialise, no clusters, and all elements belong to no clusters
  n <- nrow(sim_mat)
  clust <- list()
  spares <- rep(TRUE, n)
  
  clust_id <- 1
  ##now we keep working until spares is all false
  while(any(spares)) {
    debug_spares_loop <- sum(spares)
    debug_spares_expected <- debug_spares_loop
    new_clust <- rep(FALSE, n)
    valid_seeds <- spares

    aff_local <- rep(0, n)

    ##Do the next steps until stability is reached
    repeat{
      is_change <- FALSE



      ##find initial cluster point
      if(all(!new_clust)){
        maxima <- which(spares & valid_seeds)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        next_seed <- new_elem
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        is_change <- TRUE
        debug_spares_expected <- debug_spares_expected - 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "first"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##addition stage
      ##affinity seems to be used to bias group size
      ##Strong affinity is needed to join a large group

      ##Short-circuit evaluation, won't try to find max of an empty set
      ##trialling a slightly different approach.
      ##Find any high affinity elements.
      ##high affinity, sum of affinities into cluster (which is aff_local)
      ##exceeds thres * size of cluster
      ##So now it is a rolling mean. At each iteration, you add to the sum,
      ## then only divide by n when you need to test.
      while(any(spares) && (max(aff_local[spares]/sum(new_clust)) >= aff_thres)){
        maxima <- which(aff_local == max(aff_local[spares]) & spares)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        debug_spares_expected <- debug_spares_expected - 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "add"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##Removal stage
      ##The `-1` in the denominator excludes the self affinity,
      ##a cluster with 3 elements only has 2 affinity scores to average per element.
      ##while elements from another cluster would have 3 affinity scores to average when testing affinity to this cluster.
      while(sum(new_clust) > 1 && (min(aff_local[new_clust]/(sum(new_clust)-1)) < aff_thres)){
        minima <- which(aff_local == min(aff_local[new_clust]) & new_clust)
        new_elem <- minima[sample.int(length(minima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- FALSE
        spares[new_elem] <- TRUE
        debug_spares_expected <- debug_spares_expected + 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "remove"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        aff_local <- aff_local - aff_func(new_clust | spares, new_elem, sim_mat)
      }

      if(sum(new_clust) == 1){
        ##cluster has only one element
        if(next_seed == which(new_clust)){
          ##Since seed failed to grow beyond one element, create an outlier cluster
         break
        }
        ##If, somehow, we end up with a different final element to the seed, restart with
        ## different seed
        valid_seeds[next_seed] <- FALSE
        #message("seeds left: [", sum(valid_seeds), "]")
        if(all(!valid_seeds)){
          ##no more valid seeds exist, all have been tried, create a leftovers cluster
          new_clust <- spares
          debug_spares_expected <- debug_spares_expected - sum(new_clust)
          spares <- rep(FALSE, n)
          break
        }
      }

      if(!is_change){
        break
      }
    }
    #message("cluster assigned of size [", sum(new_clust), "]")
    #message("[", sum(spares), "] sites left to assign.")
    assertthat::assert_that(debug_spares_expected == (debug_spares_loop - sum(new_clust)))
    clust[[clust_id]] <- which(new_clust)
    clust_id <- clust_id + 1


  }

  return(clust)
}
