#' Modet-robust standardization in CRT Point Estimate
#'
#' This function calculates a model-robust point estimate for a clustered randomized trial (CRT).
#'
#' @param data A data frame where categorical variables should already be converted to dummy variables.
#' @param cluster A string representing the column name of the cluster ID in the data frame.
#' @param formula A formula for the outcome mean model, including covariates.
#' @param trt A string representing the column name of the treatment assignment per cluster.
#' @param family The link function for the outcome. Can be one of the following:
#'     - `gaussian(link = "identity")`: for continuous outcomes. Default is gaussian("identity")
#'     - `binomial(link = "logit")`: for binary outcomes.
#'     - `poisson(link = "log")`: for count outcomes.
#'     - `gaussian(link = "logit")`: for binary outcomes with logit link to model the genealized linear model.
#' @param corstr A string specifying the correlation structure for GEE models (e.g., "exchangeable", "independence").
#' @param method A string specifying the outcome mean model. Possible values are:
#'     - 'GLM': Generalized linear model on cluster-level means (continous/binary outcome).
#'     - 'LMM': linear mixed model on individual-level observations (continuous outcome).
#'     - 'GEE': marginal models fitted by generalized estimating equations.
#'     - 'GLMM': generalized linear mixed model.
#' @param trtprob A vector of treatment probabilities per cluster (for each individual), conditional on covariates. Default is rep(0.5,nrow(data))
#' @param scale A string specifying the risk measure of interest. Can be 'RD' (risk difference), 'RR' (relative risk), or 'OR' (odds ratio).
#' @return A list with the following components:
#'   - `data1`: A data frame containing all individual-level observations.
#'   - `data_clus`: A data frame contaning all cluster-level summaries.
#'   - `c(cate,iate,test_NICS)`: A vector containing: (i) cate: point estimate for cluster-average treatment effect;
#'                               (ii) iate: point estimate for individual-average treatment effect; (iii) test_NICS: value of test statistics for non-informative cluster sizes.

MRStdCRT_point <- function(formula, data, cluster, trt, trtprob,
                             family=gaussian(link="identity"),
                             corstr, method="GLM", scale){
  ################################################################
  #                                                              #
  #   Input:                                                     #
  #   formula: formula of outcome mean model.                    #
  #   data: data frame, where the categorical variables          #
  #       should be already converted to dummy variables.        #
  #   cluster: a character of the variable name of the           #
  #            cluster id.                                       #
  #   trt: a character of the variable name of the treatment     #
  #        assignment per cluster.                               #
  #   trtprob: a vector of treatment probabilities per cluster      #
  #         conditional on covariates.                           #
  #   method: specifications of outcome mean models;             #
  #     potential values are:                                    #
  #       (1) 'GLM' (linear model on cluster-level means),       #
  #       (2) 'LMM' (linear mixed model on individual-level      #
  #                 observations),                               #
  #       (3) 'GLMM' (marginal models fitted by generalized      #
  #                   linear mixed model),                       #
  #       (4) 'GEE' (generalized estimating equations).          #
  #   family: the link function for outcome                      #
  #           gaussian(link = "identity")                        #
  #           binomial(link = "logit")                           #
  #           poisson(link = "log")                              #
  #   corstr: correlation structure for GEE model                #
  #           "exchangeable","independence", etc.                #
  #   scale: risk differences ('RD'), relative risks ('RR'),     #
  #          and odds ratios ('OR').                             #
  #                                                              #
  ################################################################

  ## This file contains functions calculating point estimates, associated 95% CIs,
  ## and p-values for testing the non-informative cluster sizes for the methods
  ## described in the manuscript titled ‘Model-Robust Standardization in Cluster-Randomized Trials.’

  # Function calculating point estimates

  # Example:
  # formula = Y~X1+X2+cluster(N+H1+H2)


  # Validate inputs

  tryCatch({
    stopifnot(is.data.frame(data))
  }, error = function(e) {
    stop("Error: The provided object is not of class 'data.frame'.")
  })

  tryCatch({
    stopifnot(inherits(family, "family"))
  }, error = function(e) {
    stop("Error: The provided family is not a 'family' object.")
  })

  stopifnot(is.character(cluster), is.character(trt))

  if (method == "GLMM") {
    fam <- family$family
    lnk <- family$link

    is_allowed <-
      (fam == "gaussian"  && lnk == "identity") ||
      (fam == "binomial"  && lnk == "logit")   ||
      (fam == "poisson"   && lnk == "log")

    if (!is_allowed) {
      stop(
        "`family` for method = \"GLMM\" must be one of:\n",
        "  - gaussian(link = \"identity\")\n",
        "  - binomial(link = \"logit\")\n",
        "  - poisson(link = \"log\")\n",
        "You supplied: ", fam, "(link = \"", lnk, "\")",
        call. = FALSE
      )
    }
  }


  # Extract all variable names from the formula
  # outcome
  outcome <- all.vars(formula[[2]])
  rhs_terms <- terms(formula)

  # Extract individual-level covariates (before 'cluster')
  adj_cov <- attr(rhs_terms, "term.labels")

  # following columns need to be checked
  columns_to_check <- c(outcome, adj_cov, cluster, trt)

  # Check if all columns exist in the data
  missing_columns <- setdiff(columns_to_check, colnames(data))

  # If any columns are missing, stop and print an error message
  if (length(missing_columns) > 0) {
    stop(paste("Error: The column(s)", paste(missing_columns, collapse = ", "), "do not exist in the data."))
  }


  # Check if the treatment is assigned at the cluster-level
  non_unique_clusters <- tapply(data[[trt]], data[[cluster]], function(x) length(unique(x)) > 1)

  # Find the cluster IDs with mixed values
  mixed_clusters <- names(non_unique_clusters)[non_unique_clusters]

  # If any clusters have mixed values for 'A', print a warning with the cluster IDs
  if (length(mixed_clusters) > 0) {
    warning(paste("Warning: The following clusters have a mixture of 1 and 0", paste(mixed_clusters, collapse = ", ")))
  }
  # Create new data set for analysis
  data1 <- data %>% select(all_of(c(cluster, trt, outcome, adj_cov))) %>%
    rename(cluster = !!cluster, A = !!trt, Y = !!outcome) %>%
    arrange(cluster)
  data1$prob <- trtprob
  if (length(adj_cov) > 0) {
    adj_cov_names <- paste0("X", seq_along(adj_cov))
    names(data1)[4:(3 + length(adj_cov))] <- adj_cov_names
  }else{
    adj_cov_names <- NULL
  }

  # Calculate cluster sizes if not exists in the cluster level covariates not not exists
  data1 <- data1 %>% group_by(cluster) %>% mutate(N = n()) %>% ungroup()


  # Generate a data frame of cluster means
  data_clus <- as.data.frame(data1 %>% group_by(cluster) %>% summarise(across(everything(), mean)))



  cluster_vars <- adj_cov_names[sapply(adj_cov_names, function(nm) {
    # for each variable, check that within every cluster there's only 1 distinct (non-NA) value
    all(tapply(data1[[nm]], data1$cluster, function(v) length(unique(v[!is.na(v)])) == 1))
  })]


  ind_cov_names <- adj_cov_names[which(!(adj_cov_names %in% cluster_vars))]
  # Generate a data frame of cluster means
  data1 <- data1 %>%
    group_by(cluster) %>%
    mutate(across(all_of(ind_cov_names), ~ mean(.x), .names = "{.col}b")) %>%
    mutate(across(all_of(ind_cov_names), ~ .x - get(paste0(cur_column(), "b")), .names = "{.col}c")) %>%
    ungroup()


  # Build the formula for the individual-level model
  formulai <- paste("Y ~ A +",
                    if (length(ind_cov_names) > 0) paste(paste0(ind_cov_names, "c"), collapse = " + ") else "",
                    if (length(ind_cov_names) > 0) paste(" + ", paste(paste0(ind_cov_names, "b"), collapse = " + "), collapse = " + ") else paste(paste0(adj_cov_names), collapse = " + "),
                    if (length(cluster_vars) > 0 && (length(ind_cov_names) > 0 || length(cluster_vars) > 0)) paste(" + ", paste(cluster_vars, collapse = " + ")) else paste(cluster_vars, collapse = " + ")

  )

  # Clean up extra "+" signs
  formulai <- gsub("\\s+\\+\\s+$", "", formulai)  # Remove trailing "+" if present
  formulai <- gsub("^\\s+", "", formulai)  # Remove leading spaces

  # Build the formula for the cluster-level model
  formulac <- paste("Y ~ A +",
                    if (length(ind_cov_names) > 0) paste(ind_cov_names, collapse = " + ") else "",
                    if (length(cluster_vars) > 0 && length(ind_cov_names) > 0) paste(" + ", paste(cluster_vars, collapse = " + ")) else paste(cluster_vars, collapse = " + ")
  )

  # Clean up extra "+" signs
  formulac <- gsub("\\s+\\+\\s+$", "", formulac)  # Remove trailing "+" if present
  formulac <- gsub("^\\s+", "", formulac)  # Remove leading spaces



  ## check the method is consistent with family
  #if (method %in% c("LM", "LMM")) {
  #  if (family$family != "gaussian" && family$link != "identity") {
  #    stop("Use 'LM' or 'LMM' with 'identity' link for continuous outcomes")
  #  }
  #}

  ## Fit model
  model <- switch(method,
                  "GLM" = try(glm(formulac, data = data_clus,family=family), silent = T),
                  "LMM" = try(lme(as.formula(formulai), random = ~ 1 | cluster, data = data1), silent = T),
                  "GEE" = try(geeglm(as.formula(formulai), id = cluster, data = data1, corstr = corstr,family=family), silent=T),
                  "GLMM" =  try(glmer(paste(formulai, "+ (1 | cluster)"), data = data1, family = family),silent = T),
                  stop("Invalid method specified.")
  )


  if (inherits(model, "try-error")) {
    warning(
      "All-0 or all-1 clusters yield infinite logit with gaussian(link='logit')."
    )
    warning(
      "Proceed with non-parametric estimators"
    )
    eta <- data_clus %>%
      select(cluster) %>%
      mutate(eta1 = 0, eta0 = 0) %>%
      group_by(cluster) %>%
      as.data.frame()


  }else{
    if (method == "GLM"){
      # eta: vector containing the predicted outcome in two arms
      eta <- data_clus %>% mutate(eta1= predict(model, newdata = mutate(data_clus, A = 1), type = "response")) %>%
        mutate(eta0 = predict(model, newdata = mutate(data_clus, A = 0), type = "response")) %>%
        dplyr::select(cluster, eta1, eta0) %>%
        group_by(cluster) %>%
        as.data.frame()
    } else if (method == "LMM") {
      if (length(adj_cov_names) == 0 ) {
        eta <- data_clus %>%
          mutate(eta1 = as.matrix(cbind(rep(1,nrow(data_clus)),rep(1,nrow(data_clus)))) %*%
                   as.vector(fixef(model)[c("(Intercept)", "A")])) %>%
          mutate(eta0 = as.matrix(cbind(rep(1,nrow(data_clus)),rep(0,nrow(data_clus))))  %*%
                   as.vector(fixef(model)[c("(Intercept)", "A")])) %>%
          as.data.frame()

      }else{
        eta <- data_clus %>% mutate(eta1= as.matrix(cbind(rep(1,nrow(data_clus)),rep(1,nrow(data_clus)),data_clus[,c(grep("^X\\d$", names(data_clus)), grep("^H", names(data_clus)))]))%*%as.vector(fixef(model)[-c(grep("^X\\d+c$", names(fixef(model))))])) %>%
          mutate(eta0 = as.matrix(cbind(rep(1,nrow(data_clus)),rep(0,nrow(data_clus)),data_clus[,c(grep("^X\\d$", names(data_clus)), grep("^H", names(data_clus)))]))%*%as.vector(fixef(model)[-c(grep("^X\\d+c$", names(fixef(model))))])) %>%
          as.data.frame()
      }

    } else if (method == "GEE") {
      eta <- data1 %>% mutate(eta1= predict(model, newdata = mutate(data1, A = 1), type = "response")) %>%
        mutate(eta0 = predict(model, newdata = mutate(data1, A = 0), type = "response")) %>%
        dplyr::select(cluster, eta1, eta0) %>%
        group_by(cluster) %>%
        summarise_all(mean) %>%
        as.data.frame()
    } else if (method == "GLMM") {
      pred1 <- as.matrix(cbind(rep(1,nrow(data1)),rep(1,nrow(data1)),data1[,c(grep("^X\\d+c$", names(data1)), grep("^X\\d+b$", names(data1)), grep("^H\\d+$", names(data1)))]))%*%as.vector(fixef(model))
      pred0 <- as.matrix(cbind(rep(1,nrow(data1)),rep(0,nrow(data1)),data1[,c(grep("^X\\d+c$", names(data1)), grep("^X\\d+b$", names(data1)), grep("^H\\d+$", names(data1)))]))%*%as.vector(fixef(model))
      if(family$link=="logit"){
        # Pi: mathematical constant
        Pi <- 3.141592653589793
        eta <- data1 %>%
          mutate(
            eta1 = exp(pred1 / sqrt(3 * as.numeric(VarCorr(model)) / Pi^2 + 1)) / (1 + exp(pred1 / sqrt(3 * as.numeric(VarCorr(model)) / Pi^2 + 1))),
            eta0 = exp(pred0 / sqrt(3 * as.numeric(VarCorr(model)) / Pi^2 + 1)) / (1 + exp(pred0 / sqrt(3 * as.numeric(VarCorr(model)) / Pi^2 + 1)))
          ) %>%
          group_by(cluster) %>%
          summarise(eta1 = mean(eta1), eta0 = mean(eta0)) %>%
          as.data.frame()

      }else if(family$link=="log"){

        eta <- data1 %>%
          mutate(
            eta1 = exp(pred1 + as.numeric(VarCorr(model)) / 2),
            eta0 = exp(pred0 + as.numeric(VarCorr(model)) / 2)
          ) %>%
          group_by(cluster) %>%
          summarise(eta1 = mean(eta1), eta0 = mean(eta0)) %>%
          as.data.frame()

      }

    }

  }



  #####################Point estimates using our proposed methods###############
  mu_C1 <- data_clus$A / data_clus$prob * (data_clus$Y - eta$eta1) + eta$eta1
  mu_C0 <- (1-data_clus$A) / (1-data_clus$prob) * (data_clus$Y - eta$eta0) + eta$eta0
  cate <- switch(scale,
                 "RD" = mean(mu_C1) - mean(mu_C0),
                 "RR" = mean(mu_C1)/mean(mu_C0),
                 "OR" =  mean(mu_C1)/(1-mean(mu_C1))/mean(mu_C0)*(1-mean(mu_C0)),
                 stop("Invalid scales specified."))
  mu_I1 <- mean(data_clus$N * mu_C1)/mean(data_clus$N)
  mu_I0 <- mean(data_clus$N * mu_C0)/mean(data_clus$N)
  iate <- switch(scale,
                 "RD" = mean(mu_I1) - mean(mu_I0),
                 "RR" = mean(mu_I1)/mean(mu_I0),
                 "OR" =  mean(mu_I1)/(1-mean(mu_I1))/mean(mu_I0)*(1-mean(mu_I0)),
                 stop("Invalid scales specified."))

  test_NICS <-  cate - iate
  return(list(data1,
              data_clus,
              c(cate,iate,test_NICS)))
}







