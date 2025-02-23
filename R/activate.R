#' Calculate Counts Per Million (CPM)
#'
#' This function calculates the counts per million (CPM) for a given count matrix.
#' It normalizes the count data by scaling each column (sample) based on the total count sum for that column.
#'
#' @param expr A matrix or data frame containing raw count data, where rows represent genes and columns represent samples.
#'
#' @return A matrix with CPM values, where each column is normalized to a sum of 1,000,000.
#' @export
#'
#' @examples
#'
#' expr <- cpm(expr)
cpm <- function(expr) {
  apply(expr, 2, function(x) { x / sum(x) * 1000000 })
}

#' Normalize Counts Matrix to CPM and Apply Log Normalization (Optional)
#'
#' This function normalizes a count matrix to CPM, applies log normalization if specified,
#' and filters out rows with low mean CPM values across samples if the `filter` option is TRUE.
#'
#' @param counts.matrix A matrix or data frame containing raw count data, where rows represent genes and columns represent samples.
#' @param filter A logical value (default = TRUE). If TRUE, filters out genes with an average CPM less than 0.25 across samples.
#' @param log_normalize A logical value (default = TRUE). If TRUE, applies log transformation to the normalized CPM values (with a small offset to avoid log(0)).
#'
#' @return A matrix of CPM values after normalization and (optionally) log transformation and filtering.
#' @export
#'
#' @examples
#' expr <- cpm_normalized(expr)
cpm_normalized <- function(counts.matrix, filter = TRUE, log_normalize = TRUE) {
  loggeomeans <- rowMeans(log(counts.matrix))

  # Calculate normalization scaling factor for each sample
  sf = apply(counts.matrix, 2, function(cnts) {
    exp(stats::median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
  })

  # Normalize counts matrix based on the scaling factor
  normalized.matrix = t(t(counts.matrix) / sf)

  # Convert normalized counts to CPM
  normalized.cpm = apply(normalized.matrix, 2, function(x) { x / sum(x) * 1000000 })

  # Apply filtering if enabled
  if (filter) {
    normalized.cpm <- normalized.cpm[rowMeans(normalized.cpm) > 0.25, ]
  }

  # Apply log normalization if enabled
  if (log_normalize) {
    normalized.cpm = log(normalized.cpm + 0.001)
  }

  return(normalized.cpm)
}

#' Gene Expression Activity Estimation Using Model Matrix
#'
#' This function calculates gene activity based on a model matrix (such as pathway or gene signature models).
#' It multiplies the gene expression matrix with the model matrix for a set of common genes.
#'
#' @param expr A matrix or data frame containing gene expression data, where rows represent genes and columns represent samples.
#' @param model A matrix or data frame containing the model for gene activity (e.g., gene signature or pathway information).
#' @param scale A logical value (default = TRUE). If TRUE, scales the resulting activity scores for each sample.
#'
#' @return A matrix of gene activity scores for each sample based on the provided model matrix.
#' @export
#'
#' @examples
#' ac <- activity(expr, model, scale = FALSE)
activity <- function(expr, model, scale = TRUE) {
  # Identify common genes between the expression data and the model
  common_genes <- intersect(rownames(expr), rownames(model))

  # Calculate gene activity by multiplying the expression matrix with the model matrix
  result <- t(expr[common_genes, , drop = FALSE]) %*% as.matrix(model[common_genes, , drop = FALSE])

  # Scale the result if the scale option is enabled
  if (nrow(result) > 1 & scale == TRUE) {
    rn <- rownames(result)
    result <- apply(result, 2, scale)
    rownames(result) <- rn
  }

  return(t(result))
}
