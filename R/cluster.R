#' compute.l2fc
#' calculate l2fc between clusters
#'
#' @param data dataset (transposed)
#' @param labels cluster label values
#'
#' @return vector
#' @export
#'
#' @examples
#' Given data (a dataframe with patient samples in columns and genes in rows)
#' and a set of labels (a list of phenotypes for samples, length = ncol(data)) you can use
#' the following command to output a list of log2foldchange (length = nrow(data)) 
#' for each gene calculated between distinct phenotypes.
#' 
#' cluster_lfc <- compute.l2fc(data, labels)
compute.l2fc <- function(data, labels){
  # genes have to be in columns of data

  if (nrow(data) != length(labels)){
    stop("Rows in data not equal to length of Labels")
  }
  else if (!is.factor(labels)){
    stop("Labels is not factor")
  }
  else if (length(levels(labels)) > 2){
    stop("Factors Labels should be length 2")
  }

  ls <- levels(labels)
  l2fc <- log2(apply(data[labels == ls[2], , drop=F], 2, mean)) - log2(apply(data[labels == ls[1], , drop=F], 2, mean))

  return(l2fc)
}

#' cluster.l2fc
#' cluster data based on l2fc
#'
#' @param mat dataset
#' @param data_l2fc known l2fc
#'
#' @return list of matrices
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @export
#'
#' @examples
#' Given mat (a matrix of gene transcript counts)
#' and data_l2fc (a reference log2foldchange for each gene with gene names, 
#' calculated from the TCGA dataset using the DESeq2 library in our example),
#' this function uses hclust to cluster the new matrix counts into labels and 
#' calculates the l2fc using these new labels on mat. The function returns the labels and
#' new l2fc in a list:
#' 
#' predicted_cluster <- cluster.l2fc(mat, data_l2fc)
cluster.l2fc <- function(mat, data_l2fc){
  commongene <- intersect(rownames(mat), names(data_l2fc))
  mat <- t(mat)[, commongene]
  data_l2fc <- data_l2fc[commongene]

  # dist_mat <- dist(mat)
  dist_mat <- as.dist(1-cor(t(mat)))

  #cluster samples into primitive and committed
  hc <- hclust(dist_mat, "complete")
  model <- cutree(hc, k = 2)
  labels_cluster <- as.factor(model)

  cat(sprintf("Rows of Matrix %d\n", nrow(mat)))
  cat(sprintf("Cols of Matrix %d\n", ncol(mat)))
  cat(sprintf("Length of Cluster Labels %d\n", length(labels_cluster)))
  cat(sprintf("Length of DeSeq L2fc %d\n", length(data_l2fc)))

  # log2foldchange <- data_l2fc[names(data_l2fc) %in% colnames(mat)]
  log2foldchange_cluster <- compute.l2fc(mat, labels_cluster)

  return(list(cluster = labels_cluster, lfc = log2foldchange_cluster))
}

