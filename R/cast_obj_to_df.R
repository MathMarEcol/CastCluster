#' Convert CAST object to a data.frame mapping elements to clusters
#'
#' CAST objects are implemented as a list of clusters, with each
#' cluster storing the set of elements making up the cluster.
#' Plotting usually needs a data.frame. This helper function
#' produces a ggplot ready data.frame.
#' @param cast_ob A CAST object returned by a CAST method
#' @return data.frame with $elem and $clust cols. Elem is a row id of the data passed into the CAST method.
#'
#' @export
cast_obj_to_df <- function(cast_ob){
  clust_ind <- lapply(seq_along(cast_ob), function(clust, cast_ob){
    if(length(cast_ob[[clust]]) > 0) {
      return(data.frame(elem = cast_ob[[clust]], clust = clust))
    } else {
      return(NULL)
    }
  }, cast_ob = cast_ob)
  clust_ind <- do.call(rbind, clust_ind)
  clust_ind <- clust_ind[order(clust_ind$elem), ]
  return(clust_ind)
}
