#'Assign site to cluster
#'takes an affinity matrix,
#' rows are sites to assign to clusters,
#' and cols are sites that make up the cast_obj
#' so the dimensions of the matrix
#' are n by length(cast_obj[[1]])
#' the rows are not expected to exist in
#' cast_obj
#' If aff_thres is provided, set max to NA if all
#' clusters fall below aff_thres. Useful for
#' identifying new sites that are far from any
#' existing clusters, and may represent
#' novel environmental conditions.

#'find which clusters a site could belong, given aff_thres
#'returns a data.frame with site, cluster, affinity, and possibly probability
#' @export
predict_clust <- function(cast_obj,
                          new_sim_mat,
													aff_thres = 0,
                          type = c("max", "raw", "probability")
                          ){
    type <- match.arg(type)
    raw <- vapply(seq_along(cast_obj), function(clust, cast_obj, new_sim_mat){
				if(length(cast_obj[[clust]]) > 1){
						return(.rowMeans(new_sim_mat[, cast_obj[[clust]]], nrow(new_sim_mat), length(cast_obj[[clust]])))
				} else {
						return(new_sim_mat[, cast_obj[[clust]]])
				}
		}, numeric(nrow(new_sim_mat)), cast_obj = cast_obj, new_sim_mat = new_sim_mat)

    out <- switch(type,
           "max" = {
               ## For each x_row, find clust associated with max aff
               apply(raw, 1,
										 function(xr, cast_obj) {
												 if (all(xr < aff_thres)){
														 return(NA)
												 } else {
														 return(which.max(xr))
												 }
										 }, cast_obj = cast_obj)
           },
           "raw" = raw,
           "probability" = {
               ## Replace aff with probability of membership
               raw / rowSums(raw)
           })
		return(out)
}

#' Calculate membership matix from cast object
#'
#' Helper function for hubert_gamma(),
#'
#' Produces a matrix where element (i,j)
#' is 1 if site i and site j are in the same
#' cluster, otherwise 0.
#'
#' is_long toggles between long form and wide form
#' matrix for return values.
#' For very large areas, a full nxn matrix may
#' not fit it memory. However, the matrix becomes
#' sparse, as k increases, so long form should always work.
#'
#' @export
membership_mat <- function(cast_ob, is_long = FALSE){
  n_sites <- do.call(sum, lapply(cast_ob, length))

  mat_long  <-do.call(rbind, lapply(cast_ob, function(clust){
    mat_long_clust <- expand.grid(clust, clust)
    return(mat_long_clust)
  })
  )
  names(mat_long) <- c("i", "j")

  if(is_long){
    return(data.frame(mat_long, clust = 1))
  } else {
    ##I need side effects to avoid creating k large matricies
    mat_wide <- matrix(0, nrow = n_sites, ncol = n_sites)
    for (clust in cast_ob) {
      mat_wide[clust, clust] <- 1
    }
    return(mat_wide)
  }
}

#' Affinity of element `x` to elements `u`
aff_func <- function(x, u, sim_mat){
  ##Takes the vector x of rows and the row u of sim_match
  ##and returns the affinity from u to x
  assertthat::assert_that(length(x) == nrow(sim_mat))
  assertthat::assert_that(is.logical(x))
  assertthat::assert_that(is.numeric(u))

  ##Simple, this is just elements of sim_mat
  return(sim_mat[, u] * x)
}

#' Affinity of each element to each cluster
aff_clust_sum <- function(sim_mat, cast_ob){
  diag(sim_mat) <- 0
  do.call(cbind, lapply(seq_along(cast_ob), 
                        function(clust, sim_mat, cast_ob){

                          ##get the affinity to each cluster
                          if(length(cast_ob[[clust]]) > 0){
                            as_cols <- sim_mat[cast_ob[[clust]], , drop = FALSE]
                            n <- dim(as_cols)[1]
                            dn <- dim(as_cols)[2]
                            return(.Internal(colSums(as_cols, n, dn, TRUE)))
                          } else {
                            ##empty clusters have 0 affinity
                            return(rep(0, nrow(sim_mat)))
                          }
                        }, cast_ob = cast_ob, sim_mat = sim_mat))
}

#' Internal function for aff_clust_mean,
#'
#' This function may be called directly if you
#' can ensure that the diagonals of sim_mat
#' are set to NA  `diag(sim_mat) <- NA`
#' to speed up computation.
.aff_clust_mean <- function(sim_mat, cast_ob, aff_thres){
  aff_means <- do.call(cbind, lapply(seq_along(cast_ob),
                        function(clust, sim_mat, cast_ob){

                          ##get the affinity to each cluster
                          if(length(cast_ob[[clust]]) > 0){
                            as_cols <- sim_mat[cast_ob[[clust]], , drop = FALSE]
                            n <- dim(as_cols)[1]
                            dn <- dim(as_cols)[2]
                            return(.Internal(colMeans(as_cols, n, dn, TRUE)))
                          } else {
                            ##empty clusters have 0 affinity
                            return(rep(0, nrow(sim_mat)))
                          }
                        }, cast_ob = cast_ob, sim_mat = sim_mat))
  aff_means[is.nan(aff_means)] <- aff_thres
  return(aff_means)
  }

#' Find mean affinity from each element to each cluster
aff_clust_mean <- function(sim_mat, cast_ob, aff_thres){
  if(!all(is.na(diag(sim_mat)))) {
    diag(sim_mat) <- NA
   }
  castcluster:::.aff_clust_means(sim_mat, cast_ob, aff_thres)
}

#' Not used, behaves like aff_clust_mean
aff_clust_all <- function(sim_mat, cast_ob){
  do.call(cbind, lapply(seq_along(cast_ob),
                        function(clust, sim_mat, cast_ob){

                          ##get the affinity to each cluster
                          if(length(cast_ob[[clust]]) > 0){
                            return(rowMeans(sim_mat[, cast_ob[[clust]], drop = FALSE]))
                          } else {
                            ##empty clusters have 0 affinity
                            return(rep(0, nrow(sim_mat)))
                          }
                        }, cast_ob = cast_ob, sim_mat = sim_mat))
}

#' Across cluster affinity
#' An assymetric distance, where the distance from a to b
#' is the average affinity of elements in a to elements in b
#' returns a long form dataframe
aff_cluster_between <- function(cast_obj, sim_mat, aff_thres){
  diag(sim_mat) <- NA
  pairs <- expand.grid(seq_along(cast_obj), seq_along(cast_obj))
  affs <- apply(pairs, 1, function(p, cast_obj, sim_mat){
    elems_a <- cast_obj[[p[1]]]
    elems_b <- cast_obj[[p[2]]]
    return(mean(sim_mat[elems_a, elems_b], na.rm = TRUE))
  }, cast_obj = cast_obj, sim_mat = sim_mat)
  affs[is.nan(affs)] <- aff_thres
  return(data.frame(x = pairs[,1], y = pairs[,2], affs))
}

##within cluster affinity
##returns a list
aff_clust_inner <- function(cast_obj, sim_mat){
  lapply(seq_along(cast_obj), function(clust, cast_obj, sim_mat){
    elem_cor <- sim_mat[cast_obj[[clust]], cast_obj[[clust]] ]
    mean_aff <- mean(elem_cor)
    return(mean_aff)
  }, cast_obj = cast_obj, sim_mat = sim_mat)
}
