#' ggplot object of the similarity matrix
#'
#' Without a cast object, it will just plot them
#' similarity matrix as is, with no particular order.
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
gg_sim_mat <- function(sim_mat,
                       cast_ob = NULL,
                       sort_between = TRUE,
                       sort_within = TRUE,
                       highlight = FALSE
                       ){

  if(!is.null(cast_ob)) {
    ##order by cast object

    if(sort_within){
      cast_ob <- lapply(cast_ob, function(clust, sim_mat){
        ##sort by strength of affinity to cluster
        ord <- order(rowMeans(sim_mat[clust, clust, drop = FALSE]))
        return(clust[ord])
      },  sim_mat = sim_mat)
    }
    if(sort_between) {
      aff_btw <- aff_cluster_between(sim_mat = sim_mat, cast_obj = cast_ob)

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
      rects <- do.call("rbind", lapply(seq(1, length(cast_ob)), function(i, cast_ob_i){
        if(i > 1){
          start <- do.call(sum, lapply(1:(i-1), function(j, cast_ob_j){
            length(cast_ob_j[[j]])
          }, cast_ob_j = cast_ob_i)) + 1
        } else{
          start <- 1
        }
        end <- start + length(cast_ob_i[[i]])
        return(data.frame(xmin = start, xmax = end))
      }, cast_ob_i = cast_ob[clust_reorder]))
    }

  } else {
    site_reorder <-  seq.int(1, nrow(sim_mat))
  }
  sim_mat_reorder <- sim_mat[site_reorder, site_reorder]
  sim_mat_long <- data.frame(expand.grid(1:nrow(sim_mat), 1:ncol(sim_mat)), as.vector(as.matrix(sim_mat_reorder)))
  names(sim_mat_long) <- c("x", "y", "p")
  p <- ggplot2::ggplot(data = sim_mat_long,
                  mapping = aes(x = x, y = y, fill = p)) +
    ggplot2::geom_raster()
  if(highlight){
    p <- p + ggplot2::annotate(geom = "rect", xmin = rects$xmin, ymin = rects$xmin, xmax = rects$xmax, ymax = rects$xmax,  colour = "red", fill = NA)
  }
  return(p)
}
