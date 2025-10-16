#' Example Dataset: Simulated CRT (continuous outcome)
#'
#' A simulated dataset for demonstrating \pkg{MRStdCRT} with a continuous outcome.
#' Treatment is assigned at the cluster level and is constant within cluster.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{cluster_id}{Cluster identifier (integer or factor), constant within cluster.}
#'   \item{A}{Cluster-level treatment assignment (0/1), constant within cluster.}
#'   \item{Y}{Continuous outcome (numeric).}
#'   \item{X1}{Continuous Individual-level covariate (numeric).}
#'   \item{X2}{Binary Individual-level covariate (numeric).}
#'   \item{N}{Cluster size recorded on each row (repeats within cluster).}
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
#' @format A data frame with the following variables:
#' \describe{
#'   \item{cluster_id}{Cluster identifier (integer or factor), constant within cluster.}
#'   \item{A}{Cluster-level treatment assignment (0/1), constant within cluster.}
#'   \item{Y}{Binary outcome (0/1).}
#'   \item{X1}{Continuous Individual-level covariate (numeric).}
#'   \item{X2}{Binary Individual-level covariate (numeric).}
#'   \item{N}{Cluster size recorded on each row (repeats within cluster).}
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
