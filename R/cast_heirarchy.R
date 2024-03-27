# Copyright 2017-2024 Philip Dyer
# SPDX-License-Identifier: GPL-3.0-only
#' heirarchical grouping of cast objects
#' the aff_thres_vec will be sorted
#' this function is recursive
#' @export
cast_h <- function(sim_mat, aff_thres_vec, stabilize = FALSE, max_stab_iter = nrow(sim_mat) * 4){
  if(nrow(sim_mat == 1)) {
    return(list(list(clust = list(c(1)), affs = sim_mat)))
  }
  aff_thres_vec <- sort(aff_thres_vec, decreasing = TRUE)
  next_level <- cast_alg(sim_mat = sim_mat, aff_thres = aff_thres_vec[1])
  if(stabilize){
    next_level <- cast_stabilize(sim_mat = sim_mat, cast_obj = next_level,
                                 max_iter = max_stab_iter,
                                 aff_thres = aff_thres_vec[1])
  }
  if(length(aff_thres_vec) > 1){
    tmp <- castcluster:::aff_cluster_between(cast_obj = next_level, sim_mat = sim_mat, aff_thres = aff_thres_vec[1])
     new_sim_mat <- matrix(tmp$affs, nrow = length(next_level), ncol = length(next_level))
    return(c(cast_h(sim_mat = new_sim_mat,
                    aff_thres_vec = aff_thres_vec[-1],
                    stabilize = stabilize,
                    max_stab_iter = max_stab_iter),
             list(list(clust = next_level, affs = new_sim_mat)) ))
  } else {
    tmp <- castcluster:::aff_cluster_between(cast_obj = next_level, sim_mat = sim_mat)
    new_sim_mat <- matrix(tmp$affs, nrow = length(next_level), ncol = length(next_level))
    return(list(list(clust = next_level, affs = new_sim_mat)))
  }
}

#' unrolls the cast clusters to the desired depth
#' also recursive
#' @export
cast_h_reorder <- function(cast_h_obj, depth = length(cast_h_obj), reordering = NULL){
  assertthat::assert_that(depth > 0)
  assertthat::assert_that(depth <= length(cast_h_obj))

  if(is.null(reordering)){
    reordering <- 1:length(cast_h_obj[[1]]$clust)
  }
  layers <- seq(1, depth )
  reordered <- do.call("c", cast_h_obj[[1]]$clust[reordering])
  if(depth > 1){
    return(cast_h_reorder(cast_h_obj[-1], depth = depth - 1, reordering = reordered))
  } else {
    return(reordered)
  }
}
