#' ggplot object of the similarity matrix
#'
#' Without a cast object, it will just plot them
#' similarity matrix as is, with the order shuffled randomly.
#'
#' With a cast object, the similarity matrix will
#' be reordered to put sites in the same cluster
#' together.
#'
#' If secondary_sort is true, then the clusters will
#' be further reordered to put similar clusters together,
#' using a dendrogram sort.
#'
#' if highlight = TRUE, then rectangles will be drawn
#' around each cluster.
#'
#' If the mean similarity is less than log_trans_thres, the
#' similarity scores will be coloured on a log scale to
#' allow small similarities to be visualised more easily.
#' @export
gg_sim_mat <- function(sim_mat,
                       cast_ob = NULL,
                       aff_thres = NULL,
                       sort_among_clust = TRUE,
                       sort_within_clust = TRUE,
                       highlight = FALSE,
                       legend_label = "Similarity",
                       log_trans_thres = 1e-2
                       ){

  if(!is.null(cast_ob)) {
    ##order by cast object

    if(sort_within_clust){
      cast_ob <- lapply(cast_ob, function(clust, sim_mat){
        if(length(clust) > 1){
          aff_btw_wide <- sim_mat[clust,clust]

          aff_dist <- dist(aff_btw_wide)
          aff_sort <- hclust(aff_dist)

          clust_reorder <- aff_sort$order
          ##sort by strength of affinity to cluster
          return(clust[clust_reorder])
        } else {
          return(clust)
        }
      },  sim_mat = sim_mat)
    }
    if(sort_among_clust) {
      aff_btw <- castcluster:::aff_cluster_between(sim_mat = sim_mat, cast_obj = cast_ob, aff_thres = aff_thres)

      aff_btw_wide <- matrix(aff_btw$affs, sqrt(nrow(aff_btw)), sqrt(nrow(aff_btw)))

      aff_dist <- dist(aff_btw_wide)
      aff_sort <- hclust(aff_dist)

      clust_reorder <- aff_sort$order

    } else {
      clust_reorder <- seq(1, length(cast_ob))
    }
    site_reorder <- do.call("c", cast_ob[clust_reorder])

    ##add cluster rectangles
    if(highlight) {
      offset <- 0.5
      rects <- do.call("rbind", lapply(seq(1, length(cast_ob)), function(i, cast_ob_i){
        if(i > 1){
          start <- do.call(sum, lapply(1:(i-1), function(j, cast_ob_j){
            length(cast_ob_j[[j]])
          }, cast_ob_j = cast_ob_i)) + 1
        } else{
          start <- 1
        }
        end <- start + length(cast_ob_i[[i]])
        return(data.frame(xmin = start - offset, xmax = end - offset))
      }, cast_ob_i = cast_ob[clust_reorder]))
    }

  } else {
    site_reorder <-  sample.int(nrow(sim_mat))
  }
  sim_mat_reorder <- sim_mat[site_reorder, site_reorder]
  sim_mat_long <- data.frame(expand.grid(1:nrow(sim_mat), 1:ncol(sim_mat)), as.vector(as.matrix(sim_mat_reorder)))
  names(sim_mat_long) <- c("x", "y", "p")
  p <- ggplot2::ggplot(data = sim_mat_long,
                  mapping = ggplot2::aes(x = x, y = y, fill = p)) +
    ggplot2::geom_raster() +
    ggplot2::coord_fixed() +
    ggthemes::theme_tufte() +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_continuous(position = "top") +
    ggplot2::labs(fill = legend_label, x = NULL, y = NULL) +
    if ( mean(sim_mat[lower.tri(sim_mat)]) < log_trans_thres ) {
      ggplot2::scale_fill_gradient(low = "black", high = "white", breaks = c(0,log_trans_thres,1), trans = scales::pseudo_log_trans(base = 10, sigma = log_trans_thres/(nrow(sim_mat)*(nrow(sim_mat)-1)/2) ))
    } else {
      ggplot2::scale_fill_gradient(low = "black", high = "white")
    }

  if(highlight & !is.null(cast_ob)){
    p <- p + ggplot2::annotate(geom = "rect", xmin = rects$xmin, ymin = rects$xmin, xmax = rects$xmax, ymax = rects$xmax,  colour = "red", fill = NA)
  }
  return(p)
}

#' Get sort order of clusters, so similar clusters are put together
#'
#' @export
sort_between <- function(sim_mat, cast_ob){

  aff_btw <- castcluster:::aff_cluster_between(sim_mat = sim_mat, cast_obj = cast_ob)

  aff_btw_wide <- matrix(aff_btw$affs, sqrt(nrow(aff_btw)), sqrt(nrow(aff_btw)))

  aff_dist <- dist(aff_btw_wide)
  aff_sort <- hclust(aff_dist)

  return(aff_sort$order)

}

