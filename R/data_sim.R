#' Example Dataset: Simulated CRT (continuous outcome)
#'
#' A simulated dataset for demonstrating \pkg{MRStdCRT} with a continuous outcome.
#' Treatment is assigned at the cluster level and is constant within cluster.
#'
#' @format A data frame with the following variables (10 columns):
#' \describe{
#'   \item{A}{Cluster-level treatment assignment (0/1), constant within cluster.}
#'   \item{H1}{Cluster-level covariate 1.}
#'   \item{H2}{Cluster-level covariate 2.}
#'   \item{N}{Cluster size recorded on each row (repeats within cluster).}
#'   \item{X1}{Individual-level covariate 1 (numeric).}
#'   \item{X2}{Individual-level covariate 2 (numeric or binary coded 0/1).}
#'   \item{Y}{Observed continuous outcome.}
#'   \item{Y0}{Potential outcome under control (continuous).}
#'   \item{Y1}{Potential outcome under treatment (continuous).}
#'   \item{cluster_id}{Cluster identifier (integer or factor), constant within cluster.}
#' }
#'
#' @usage data(data_sim_continuous)
#' @keywords datasets
#' @source Simulated data included with the package for examples.
#' @examples
#' data(data_sim_continuous)
#' head(data_sim_continuous)
#' table(data_sim_continuous$cluster_id)
"data_sim_continuous"


#' Example Dataset: Simulated CRT (binary outcome)
#'
#' A simulated dataset for demonstrating \pkg{MRStdCRT} with a binary outcome.
#' Treatment is assigned at the cluster level and is constant within cluster.
#'
#' @format A data frame with the following variables (10 columns):
#' \describe{
#'   \item{A}{Cluster-level treatment assignment (0/1), constant within cluster.}
#'   \item{H1}{Cluster-level covariate 1.}
#'   \item{H2}{Cluster-level covariate 2.}
#'   \item{N}{Cluster size recorded on each row (repeats within cluster).}
#'   \item{X1}{Individual-level covariate 1 (numeric).}
#'   \item{X2}{Individual-level covariate 2 (numeric or binary coded 0/1).}
#'   \item{Y}{Observed binary outcome (0/1).}
#'   \item{Y0}{Potential outcome under control (0/1).}
#'   \item{Y1}{Potential outcome under treatment (0/1).}
#'   \item{cluster_id}{Cluster identifier (integer or factor), constant within cluster.}
#' }
#'
#' @usage data(data_sim_binary)
#' @keywords datasets
#' @source Simulated data included with the package for examples.
#' @examples
#' data(data_sim_binary)
#' head(data_sim_binary)
#' with(data_sim_binary, table(A, Y))
"data_sim_binary"
