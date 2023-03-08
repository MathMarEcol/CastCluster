#' Two functions that extend cast to reduce the
#' number of clusters
#'
#' cast_alg_cautious stabilizes the matrix after every cluster is added,
#' so we don't create a lot of small clusters on the edges.
#' I expect to use all the sites faster.
cast_alg_cautious <- function(sim_mat, aff_thres, max_iter = 20){
##rather than wrestle with integer sets, I will use boolean vectors
  
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
      while(any(new_clust) && (min(aff_local[new_clust]/sum(new_clust)) < aff_thres)){
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

      if(all(!new_clust)){
        ##cluster is empty
        valid_seeds[next_seed] <- FALSE
        #message("seeds left: [", sum(valid_seeds), "]")
        if(all(!valid_seeds)){
          ##no more valid seeds exist, all have been tried, create a leftovers cluster
          new_clust <- spares
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

    ##now, stabilize the clusters.
    clust[[clust_id]] <- which(new_clust)
    if(any(spares)){
      clust[[clust_id + 1]] <- which(spares)
    }
    clust <- cast_stabilize(cast_obj = clust, aff_thres = aff_thres, sim_mat = sim_mat, max_iter = max_iter)
    if(any(spares)){
      spares <- rep(FALSE, n)
      spares[clust[[clust_id +1]] ] <- TRUE
    }
    #message("[", sum(spares), "] sites left to assign after stabilize step.")
    clust_id <- clust_id + 1
  }

  return(clust)
}
