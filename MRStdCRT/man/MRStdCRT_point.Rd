% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRStdCRT_point.R
\name{MRStdCRT_point}
\alias{MRStdCRT_point}
\title{Modet-robust standardization in CRT Point Estimate}
\usage{
MRStdCRT_point(
  formula,
  data,
  clus_id,
  trt,
  prob,
  family = gaussian(link = "identity"),
  corstr,
  method = "LM",
  scale
)
}
\arguments{
\item{formula}{A formula for the outcome mean model, including covariates.}

\item{data}{A data frame where categorical variables should already be converted to dummy variables.}

\item{clus_id}{A string representing the column name of the cluster ID in the data frame.}

\item{trt}{A string representing the column name of the treatment assignment per cluster.}

\item{prob}{A vector of treatment probabilities per cluster, conditional on covariates.}

\item{family}{The link function for the outcome. Can be one of the following:
- `gaussian(link = "identity")`: for continuous outcomes. Default is gaussian("identity")
- `binomial(link = "logit")`: for binary outcomes.
- `poisson(link = "log")`: for count outcomes.}

\item{corstr}{A string specifying the correlation structure for GEE models (e.g., "exchangeable", "independence").}

\item{method}{A string specifying the outcome mean model. Possible values are:
- 'LM': linear model on cluster-level means (continuous outcome).
- 'LMM': linear mixed model on individual-level observations (continuous outcome).
- 'GEE': marginal models fitted by generalized estimating equations.
- 'GLMM': generalized linear mixed model.}

\item{scale}{A string specifying the risk measure of interest. Can be 'RD' (risk difference), 'RR' (relative risk), or 'OR' (odds ratio).}
}
\value{
A list with the following components:
  - `data1`: A data frame containing all individual-level observations.
  - `data_clus`: A data frame contaning all cluster-level summaries.
  - `c(cate,iate,test_NICS)`: A vector containing: (i) cate: point estimate for cluster-average treatment effect;
                              (ii) iate: point estimate for individual-average treatment effect; (iii) test_NICS: value of test statistics for non-informative cluster sizes.
}
\description{
This function calculates a model-robust point estimate for a clustered randomized trial (CRT).
}
