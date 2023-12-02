source("./R/cluster.R")
source("./R/dist.R")

#' fitness.cor
#' Correlation-based fitness function
#'
#' @param individual individual from GA object
#' @param penalty.function user-specified penalty
#' @param k clusters
#' @param data original dataset
#' @param known_lfc lfc for comparison
#' @param method correlation method
#' @param ... additional params
#'
#' @return numeric value
#' @importFrom stats as.dist
#' @export
#'
#' @examples
#' Given a matrix 'individual' (a self object of the current genetic algorithm matrix), 
#' an integer 'k' (number of clusters in the data), 
#' a matrix 'data' (a dataset), and a 'method' of correlation classification,
#' this function calculates outputs the correlation between the log2foldchange
#' calculated on the predicted clusters and the known log2foldchange.
#' 
#' Ex:
#' fitness <- fitness.cor(individual, k=2, data, known_lfc, method = 'pearson')
fitness.cor <- function(individual, penalty.function = NULL, k, data, known_lfc, method = c('pearson', 'spearman', 'kendall'), ...){

  dims <- length(individual)/k
  fitness.value <- NA

  m.individual <- matrix(individual, nrow = k, ncol = dims)

  method <- match.arg(method)
  which.dists <- dist.corr(data, d2=m.individual, method = method)

  which.dists <- as.factor(which.dists)

  # to avoid a convergence for configurations different from user-specified k
  if (length(unique(which.dists)) < k) {
    # maximum penalty
    fitness.value = -1

  } else {
    cluster_lfc <- compute.l2fc(data, which.dists)
    fitness.value <- abs(cor(cluster_lfc, known_lfc))
  }

  return(fitness.value)
}
