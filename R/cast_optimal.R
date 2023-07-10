#' Find optimal clustering using Hubert's Gamma Statistic
#'
#' Given a similarity matrix `sim_mat`, where 1 is identical,
#' and 0 is completely dissimilar, finds an affinity
#' threshold that maximises Hubert's \eqn{\Gamma} statistic
#' (Tseng and Kao, 2003).
#'
#' It is possible that Hubert's \eqn{\Gamma} statistic is
#' is optimised at aff_thres = mean(sim_mat) when
#' the diagonals are ignored, which would make
#' this function redundant.
#' NOTE: I did find a counter-example where the mean
#' affinity was not the peak Gamma statistic.
#'
#' @param sim_mat A square numeric matrix with values
#' in the range [0, 1]. 0 indicates maximimum
#' dissimilarity, 1 indicates perfect similarity.
#' The Diagonals are assumed to be 1, each element
#' should be identical to itself.
#' @param return_full if TRUE, return a data.frame
#' with all metadata and cast objects. If FALSE,
#' just return row with best clustering.
#' @param aff_range minimum and maximum aff_thres
#' to test. Best left as default to test the full
#' range.
#' @param npar number of values of aff_thres
#' to test each run. Must be >=3, but very large values
#' will slow down algorithm with no improvement in
#' accuracy.
#' @param min_range if aff_range is narrower than
#' `min_range`, stop and return the best aff_thres.
#' @param min_tol if the difference in Hubert's \eqn{\Gamma}
#' statistic scores between the worst and best aff_thres
#' for a run are less than `min_tol`, stop and
#' return the best aff_thres
#'
#' @export
cast_optimal <- function(sim_mat,
                         return_full = TRUE,
                         aff_range = range(sim_mat[upper.tri(sim_mat)]),
                         m = 4,
                         min_range = diff(aff_range)/100,
                         min_tol = 0.0001) {
  ##Check input is meaningful
  assertthat::assert_that(assertthat::are_equal(dim(sim_mat)[1], dim(sim_mat)[2])) ##sim_mat must be square
  assertthat::assert_that(m >= 3)
  assertthat::assert_that(all(min_range > 0,  min_range < 1))
  assertthat::assert_that(all(min(aff_range) >= 0,  max(aff_range) <= 1, length(aff_range) == 2))
  assertthat::assert_that(min_tol > 0)

  ##begin recursion
  ret <- castcluster:::cast_optimal_recurse(sim_mat = sim_mat,
                         aff_range = aff_range,
                         m = m,
                         min_range = min_range,
                         min_tol = min_tol,
                         rec_data = NULL,
                         rec_depth = 1)

  if(return_full) {
    return(ret)
  } else {
    return(ret[which.max(ret$gamma),])
  }
}

#' Recursive calculation of optimal aff_thres
#'
#' `cast_opminal` is a wrapper function that
#' hides recursive parameters from the user.
cast_optimal_recurse <- function(sim_mat,
                         aff_range,
                         m,
                         min_range,
                         min_tol,
                         rec_data,
                         rec_depth) {

  ## Check for errors in recursion logic, or bad inputs
  assertthat::assert_that(diff(aff_range) > min_range)


  ## Calculate Hubert's \Gamma statistic for each partition
  aff_thres_parts <- seq(min(aff_range), max(aff_range), length.out = m)
  gamma_score <- do.call(rbind,
                         future.apply::future_lapply(aff_thres_parts, function(aff_thres, sim_mat, rec_depth) {
                           clust_first_pass <- castcluster::cast_alg(sim_mat, aff_thres)
                           clust_stabilise <- castcluster::cast_stabilize(clust_first_pass,
                                                                          aff_thres,
                                                                          sim_mat)

                           if(length(clust_stabilise) == 1) {
                             h <- NA
                           } else {
                             mem_mat <- castcluster::membership_mat(clust_stabilise)
                             h <- castcluster::hubert_gamma(sim_mat, mem_mat, norm_z = TRUE)
                           }
                           return(data.frame(aff_thres = aff_thres, gamma = h,
                                             k = length(clust_stabilise),
                                             cast_ob = I(list(clust_stabilise)),
                                             rec_depth = rec_depth))
                         },  sim_mat = sim_mat, rec_depth = rec_depth)
  )

  ## Find the best gamma score
  max_score_ind <- which.max(gamma_score$gamma)

  new_range <- c(0,0)
  ## if(max_score_ind > 1) {
  ##   new_range[1] <- aff_thres_parts[max_score_ind - 1]
  ## } else {
  ##   new_range[1] <- aff_thres_parts[max_score_ind]
  ## }
  ## if(max_score_ind < m) {
  ##   new_range[2] <- aff_thres_parts[max_score_ind + 1]
  ## } else {
  ##   new_range[2] <- aff_thres_parts[max_score_ind]
  ## }

  ## Remove the highest and lowest scores, unless
  ## the max is on one end
  gamma_na <- is.na(gamma_score$gamma)
  if(any(gamma_na)) {
    new_range[1] <- min(aff_thres_parts[!gamma_na])
    new_range[2] <- max(aff_thres_parts[!gamma_na])
  } else {
    new_range[1] <- aff_thres_parts[2]
    new_range[2] <- aff_thres_parts[m - 1]
    if(max_score_ind == 1) {
      new_range[1] <- min(aff_range) - diff(aff_thres_parts)[1] / 2
      if(new_range[1] < 0) {
        new_range[1] <- 0
      }
    } else if(max_score_ind == m) {
      new_range[2] <- max(aff_range) + diff(aff_thres_parts)[1] / 2
      if(new_range[1] > 1) {
        new_range[1] <- 1
      }
    }
  }

  max_aff <- aff_thres_parts[max_score_ind]

  ## Check whether to keep narrowing, or return
  if(diff(range(gamma_score$gamma, na.rm = TRUE)) < min_tol | diff(new_range) < min_range) {
    return(rbind(rec_data, gamma_score))
  } else {
    return(cast_optimal_recurse(sim_mat = sim_mat,
                        aff_range = new_range,
                        m = m,
                        min_range = min_range,
                        min_tol = min_tol,
                        rec_data = rbind(rec_data, gamma_score),
                        rec_depth = rec_depth + 1
                        )
           )
  }

  }

#' Huberts Gamma statistic
#'
#' Calculates Huberts Gamma statistic
#'
#' norm_z rescales both x and y
#' by subtracting the mean and dividing by the
#' standard deviation.
#' The norm_z variant is from Tseng and Kao 2003.
#'
#' The original Huberts Gamma statistic (Hubert and Levin 1976)
#' uses x and y as is.
#'
#' Zhao et al. 2006 proposes another variant, where
#' y is defined as the distance between the cluster
#' centers of site i and site j.
#' This variant is harder to calculate and gives
#' a knee rather than a peak.
#' This variant is not implemented.
#'
#' This version requires a full matrix for both sim_mat and
#' member_mat. Using other forms is much more complicated.
#'
#' @param sim_mat a square matrix equivalent, all values in the range 0 to 1
#' @param member_mat square matrix equivalent, 0 if a pair of sites are in different
#' clusters, 1 if the pair are in the same cluster.
#' @param norm_z whether to normalize the sim_mat and member_mat into z-scores. See Tseng and Kao 2003
#'
#' @return the Hubert Gamma score
#'
#' @export
hubert_gamma <- function(sim_mat, member_mat, norm_z = TRUE){
  assertthat::assert_that(assertthat::are_equal(dim(sim_mat), dim(member_mat)))
  diag(sim_mat) <- 1
  if(norm_z){
    mean_s <- mean(sim_mat, na.rm = TRUE)
    sd_s <- sd(sim_mat, na.rm = TRUE)
    mean_m <- mean(member_mat)
    sd_m <- sd(member_mat)

    member_mat <- (member_mat - mean_m) / sd_m
    sim_mat <- (sim_mat - mean_s) / sd_s
  } else {
    mean_m <- 0
    sd_m <- 1
    mean_s <- 0
    sd_s <- 1
  }

  h_gamma <- (1/sum(upper.tri(sim_mat))) *
    sum((upper.tri(sim_mat)*sim_mat) *
          member_mat)
  return(h_gamma)
}
