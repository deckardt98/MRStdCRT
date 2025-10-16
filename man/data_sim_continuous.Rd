
\docType{data}
\name{data_sim_continuous}
\alias{data_sim_continuous}
\title{Example Dataset: Simulated CRT (continuous outcome)}
\format{
A data frame with the following variables:
\describe{
  \item{cluster_id}{Cluster identifier (integer or factor), constant within cluster.}
  \item{A}{Cluster-level treatment assignment (0/1), constant within cluster.}
  \item{Y}{Continuous outcome (numeric).}
  \item{X1}{Continuous Individual-level covariate (numeric).}
  \item{X2}{Binary Individual-level covariate (numeric).}
  \item{N}{Cluster size recorded on each row (repeats within cluster).}
}
}
\source{
Simulated data included with the package for examples.
}
\usage{
data(data_sim_continuous)
}
\description{
A simulated dataset for demonstrating \pkg{MRStdCRT} with a continuous outcome.
Treatment is assigned at the cluster level and is constant within cluster.
}
\examples{
data(data_sim_continuous)
head(data_sim_continuous)
table(data_sim_continuous$cluster_id)
}
\keyword{datasets}
