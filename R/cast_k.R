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
      ##Merge two most similar clusters
      clust_sim <- clust_sim[!(clust_sim$Var1 == clust_sim$Var2),]
      clust_merge <- clust_sim[which.max(clust_sim$affs),]
      cast_ob[[clust_merge$Var1]] <- c(cast_ob[[clust_merge$Var1]], cast_ob[[clust_merge$Var2]])
      cast_ob[[clust_merge$Var2]] <- NULL

      ##Stabilise
      cast_ob <- cast_stabilize(cast_ob, aff_thres = 0, sim_mat)
    } else {
      ##Seed the weakest object in the weakest cluster into a new cluster
      clust_sim <- clust_sim[(clust_sim$Var1 == clust_sim$Var2),]
      clust_split <- clust_sim[which.min(clust_sim$affs),1]
      affs <- aff_clust_all(sim_mat, cast_ob)
      new_seed <- which(affs[ , clust_split] ==  min(affs[cast_ob[[clust_split]],clust_split]))
      cast_ob[[clust_split]] <- cast_ob[[clust_split]][cast_ob[[clust_split]] != new_seed]
      cast_ob[[length(cast_ob) + 1]] <- c(new_seed)

      ##Stabilise
      cast_ob <- cast_stabilize(cast_ob, aff_thres = 0, sim_mat)

    }

  }

  return(cast_ob)

}
