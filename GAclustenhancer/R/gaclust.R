# libraries are used along the code by explicit refer syntax package::function()
source("./R/fitness.R")
source("./R/dist.R")

ga.clust.env <- new.env(parent = emptyenv())

setClass(Class = "ga.clust",
         slots = c(original.data = "data.frame",
                   centers = "data.frame",
                   cluster = "vector",
                   correlation = "ANY",
                   lfc = "vector",
                   call = "ANY"))

#' populate
#' Generates an initial population of candidate centers for partitions
#'
#' @param object initial gaclust object
#' @param dims number of columns of the dataset
#' @param k number of clusters
#'
#' @return numeric matrix of size popSize by dims * k
#' @importFrom stats runif
#' @export
#'
#' @example
#' Given a GA object declared with params popSize, lower, upper:
#' object@population = populate(object, dims, k)
populate <- function(object, dims, k) {
  population <- matrix(as.double(NA), nrow = object@popSize, ncol = dims * k)

  for(j in 1:length(object@lower)){
    population[,j] <- runif(object@popSize, object@lower[j], object@upper[j])
  }

  return(population)
}

#' ga.clust
#' Main function, calls GA function to cluster data using correlation-based fitness
#'
#' @param dataset dataset to cluster
#' @param crossover.rate probability crossover occurs
#' @param k number of clusters
#' @param mutation.rate probability of mutation
#' @param elitism % of population retained in next iter
#' @param pop.size initial population size
#' @param generations repetitions
#' @param seed.ga seed
#' @param method correlation method
#' @param penalty.function fitness penalty
#' @param plot.internals plot internal messages?
#' @param known_lfc external lfc exists?
#' @param ... additional params
#'
#' @return ga.clust object
#' @export
#' @importFrom stats cor
#' @importFrom utils head
#'
#' @example
#' Given a genetic dataset DATA (ie., the rows are samples and the columns are genes), and a known log to fold change of the same genes, run:
#' DATA <- t(DATA)
#' obj_gene <- ga.clust(dataset = data.frame(DATA), k = 2,
#'             plot.internals = FALSE, seed.ga = 42,
#'             pop.size = nrow(data.frame(DATA)),
#'             known_lfc = input_l2fc)

ga.clust <- function(dataset = NULL, crossover.rate = 0.9, k = 2, #scale = FALSE,
                     mutation.rate = 0.01, elitism = 0.05, pop.size = 25,
                     generations = 100, seed.ga = 42, method = "pearson",
                     penalty.function = NULL,
                     plot.internals = TRUE,
                     known_lfc = NULL, ...) {

  # --- arguments validation --- #

  if (is.null(dataset)){
    stop("'dataset' can not be NULL")
  }

  if (!inherits(dataset, 'data.frame')) {
    stop("'dataset' must be a data.frame object.")
  }

  if (is.null(k)){
    stop("'k' can not be NULL")
  }

  if (is.numeric(k)) {
    # forces k to be an integer value
    k <- as.integer(k)

    if (k < 2) {
      stop("'k' must be a positive integer value greater than one (k > 1)")
    }

  } else if (is.character(k)) {
      #or one of the methods to estimate it: 'minimal' or 'broad'."
      stop("'k' must be a positive integer value greater than one (k > 1)")
  }

  if (!is.null(penalty.function)) {
    if (!is.function(penalty.function)) {
      stop("'penalty.function' must be a valid R function.")
    }
  }

  # --- final of arguments validation --- #

  call <- match.call()

  dims <- ncol(dataset)
  elitism.rate = floor(pop.size * elitism)

  # distance matrix
  d <- dist.corr(dataset, method = "euclidean")
  d2 <- d^2

  lowers <- apply(dataset, 2, min)
  uppers <- apply(dataset, 2, max)

  lower_bound <- unlist(lapply(lowers, function (x) { rep(x, k) } ))
  upper_bound <- unlist(lapply(uppers, function (x) { rep(x, k) } ))

  set.seed(seed.ga)

  # call GA functions
  cors <- list()
  genetic <- GA::ga(type = "real-valued",
                    seed = seed.ga,
                    population = function(object) populate(object, dims, k),
                    selection = "gareal_lrSelection",
                    mutation = "gareal_nraMutation",
                    crossover = "gareal_blxCrossover",
                    popSize = pop.size,
                    elitism = elitism.rate,
                    pmutation = mutation.rate,
                    pcrossover = crossover.rate,
                    maxiter = generations,
                    fitness = function(individual, penalty.function) fitness.cor(individual, penalty.function, k, dataset, known_lfc, method),
                    lower = lower_bound,
                    upper = upper_bound,
                    parallel = FALSE,
                    monitor = F)

  num_solutions = length(genetic@solution)/(k*dims)

  if (num_solutions == 1) {
    solution <- matrix(genetic@solution, nrow = k, ncol = dims)
    print(head(solution))
  } else {
    # if there is more than a single solution (they are identical,
    # and must be close for centroids values)
    solution <- matrix(genetic@solution[1,], nrow = k, ncol = dims)
  }

  # calculates the distance between each cluster and the data and returns the min value, between 1 and 0
  which.dists <- dist.corr(dataset, d2=solution, method = "pearson")
  print(length(which.dists))

  # computes lfc and correlation for each cluster label
  log2foldchange <- compute.l2fc(dataset, as.factor(which.dists))
  corr <- cor(log2foldchange, known_lfc)

  # builds the solution object
  solution.df <- as.data.frame(solution)
  colnames(solution.df) <- colnames(dataset)
  solution.df <- solution.df[with(solution.df, order(apply(solution.df, 1, sum))), ]

  object <- methods::new("ga.clust",
                         original.data = as.data.frame(dataset),
                         centers = solution.df,
                         cluster = as.vector(which.dists),
                         correlation = abs(corr),
                         lfc = log2foldchange,
                         call = call)


  # plot the results
  if (plot.internals) {
    plot(genetic, main = "Evolution")
    lim <- 500

    if (nrow(dataset) > lim) {
      cat("\nIMPORTANT!!!\nThe dataset contains ", nrow(dataset), "rows. To improve the quality of the graph, ga.clust will generate a file 'cor_ga_clust.pdf' in working directory.\n")
      grDevices::pdf(file = 'cor_ga_clust.pdf')
    }

    if (nrow(dataset) > lim) {
      garbage <- grDevices::dev.off()
    }
  }

  # return an object of class 'ga.clust'
  return(object)

}

#' print.ga.clust
#' Print out GA Clustering results
#'
#' @param x ga.clust object
#' @param ... additional params
#'
#' @return NULL
#' @export print.ga.clust
#' @export
#'
#' @examples
#' Given an gene object output from the ga.clust function, run:
#' print.ga.clust(GACLUST_GENE_OBJECT)
print.ga.clust <- function(x, ...) {

  cat("\n Description of GA clustering class objects: :\n")

  cat("\nOriginal data (first rows):\n")
  print(head(x@original.data))
  cat("\nCluster Centers:\n")
  print(x@centers)
  cat("\nCluster partitions:\n")
  print(x@cluster)
  cat("\nlog2FoldChange:\n")
  print(x@lfc)
  cat("\nCorrelation:\n")
  print(x@correlation)

  cat("\nCall:\n")
  print(x@call)
}

setMethod("print", "ga.clust", print.ga.clust)
