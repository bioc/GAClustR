#' dist.corr
#' calculate the distance between two datasets or distance of one dataset, if d2 is null
#'
#' @param data first dataset
#' @param d2 second dataset
#' @param method method of distance calculation
#'
#' @return distance matrix
#' @importFrom stats dist
#' @export
#'
#' @examples
#' Given a dataset data (of count data) and an optional second dataset d2, 
#' as well as a method of computing distance (ex., pearson, spearman, euclidean, etc.),
#' output 1 minus the correlation between datasets/dataset and itself if "pearson", "spearman", etc. are specified, or output the distance between
#' datasets/dataset and itself if options "euclidean", "manhattan", etc. are specified.
#' 
#' Example 1
#' correlation <- dist.corr(dataset, d2=dist, method = "pearson")
#' 
#' Example 2
#' distance <- dist.corr(dataset, method = "euclidean")
dist.corr <- function(data, d2=NULL, method=c("pearson", "kendall", "spearman", "euclidean", "manhattan")){

  if (is.null(data)){
    stop("'dataset' can not be NULL")
  }

  if (!inherits(data, 'data.frame')){
    stop("'dataset' must be a data.frame object.")
  }

  method <- match.arg(method)

  if(!is.null(d2)){
    if(method %in% c("pearson", "spearman", "kendall")){
      dists <- 1 - cor(t(data), t(d2), method = method)
    }
    else{
      dists <- Rfast::dista(data, d2, type = method, square = TRUE)
    }
    return(apply(dists, 1, which.min))
  }
  else{
    if(method %in% c("pearson", "spearman", "kendall")){
      dists <- 1 - cor(t(data), method = method)
    }
    else{
      dists <- dist(data, method = method, diag = FALSE, upper = FALSE)
    }
    return(dists)
  }
}
