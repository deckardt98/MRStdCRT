#' Model-robust standardization in CRT Point Estimate
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
  
  rhs_terms   <- terms(formula)
  term_labels <- attr(rhs_terms, "term.labels")
  
  all_vars_in_formula <- setdiff(all.vars(formula), outcome)
  
  rhs_syms <- setdiff(all_vars_in_formula, trt)
  
  trt_pat <- paste0("(^|:)", trt, "(:|$)")
  int_labels <- term_labels[grepl(trt_pat, term_labels)]
  if (length(int_labels) > 0) {
    partners_raw <- unique(unlist(strsplit(int_labels, ":", fixed = TRUE)))
    
    int_partners <- setdiff(partners_raw, trt)
    
    int_partners <- intersect(int_partners, rhs_syms)
  } else {
    int_partners <- character(0)
  }
  
  # following columns need to be checked
  columns_to_check <- c(outcome, rhs_syms, cluster, trt)
  
  missing_columns <- setdiff(columns_to_check, colnames(data))
  if (length(missing_columns) > 0) {
    stop(paste("Error: The column(s)", paste(missing_columns, collapse = ", "), "do not exist in the data."))
  }
  
  data1 <- data %>%
    dplyr::select(dplyr::all_of(c(cluster, trt, outcome, rhs_syms))) %>%
    dplyr::rename(cluster = !!cluster, A = !!trt, Y = !!outcome) %>%
    dplyr::arrange(cluster)
  
  data1$prob <- trtprob
  
  if (length(rhs_syms) > 0) {
    adj_cov_names <- paste0("X", seq_along(rhs_syms))
    map_old2new <- stats::setNames(adj_cov_names, rhs_syms)
    names(data1)[4:(3 + length(rhs_syms))] <- adj_cov_names
  } else {
    adj_cov_names <- NULL
    map_old2new   <- character(0)
  }
  
  data1 <- data1 %>% dplyr::group_by(cluster) %>% dplyr::mutate(N = dplyr::n()) %>% dplyr::ungroup()
  
  data_clus <- as.data.frame(data1 %>% dplyr::group_by(cluster) %>% dplyr::summarise(dplyr::across(dplyr::everything(), mean)))
  
  cluster_vars <- adj_cov_names[sapply(adj_cov_names, function(nm) {
    all(tapply(data1[[nm]], data1$cluster, function(v) length(unique(v[!is.na(v)])) == 1))
  })]
  ind_cov_names <- setdiff(adj_cov_names, cluster_vars)
  
  if (length(ind_cov_names) > 0) {
    data1 <- data1 %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(ind_cov_names), ~ mean(.x), .names = "{.col}b")) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(ind_cov_names), ~ .x - get(paste0(cur_column(), "b")), .names = "{.col}c")) %>%
      dplyr::ungroup()
  }
  
  partners_new <- unname(map_old2new[int_partners])
  partners_new <- partners_new[!is.na(partners_new)]
  
  main_terms_i <- c(
    "A",
    if (length(ind_cov_names) > 0) paste0(ind_cov_names, "c"),
    if (length(ind_cov_names) > 0) paste0(ind_cov_names, "b"),
    if (length(cluster_vars)    > 0) cluster_vars
  )
  
  
  ## interaction terms
  ia_terms_i <- character(0)
  if (length(partners_new) > 0) {
    ia_terms_i <- c(
      if (length(intersect(partners_new, ind_cov_names)) > 0)
        c(paste0("A:", intersect(partners_new, ind_cov_names), "c"),
          paste0("A:", intersect(partners_new, ind_cov_names), "b")),
      if (length(intersect(partners_new, cluster_vars)) > 0)
        paste0("A:", intersect(partners_new, cluster_vars))
    )
  }
  formulai <- stats::as.formula(
    paste("Y ~", paste(c(main_terms_i, ia_terms_i), collapse = " + "))
  )
  
  main_terms_c <- c(
    "A",
    if (length(ind_cov_names) > 0) ind_cov_names,
    if (length(cluster_vars)  > 0) cluster_vars
  )
  
  ia_terms_c <- character(0)
  if (length(partners_new) > 0) {
    ia_terms_c <- paste0("A:", partners_new)
  }
  formulac <- stats::as.formula(
    paste("Y ~", paste(c(main_terms_c, ia_terms_c), collapse = " + "))
  )
  
  
  # Clean up extra "+" signs
  #formulac <- gsub("\\s+\\+\\s+$", "", formulac)  # Remove trailing "+" if present
  #formulac <- gsub("^\\s+", "", formulac)  # Remove leading spaces
  
  
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
      # level = 0: use fixed effects
      nd1 <- data1; nd1$A <- 1
      nd0 <- data1; nd0$A <- 0
      
      # fixed effects coefficient
      beta <- fixef(model)
      
      # design matrix
      mm1 <- model.matrix(formulai, data = nd1)
      mm0 <- model.matrix(formulai, data = nd0)
      
      p1 <- as.numeric(mm1 %*% beta)
      p0 <- as.numeric(mm0 %*% beta)
      
      eta <- dplyr::tibble(cluster = data1$cluster, eta1 = p1, eta0 = p0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(eta1 = mean(eta1), eta0 = mean(eta0), .groups = "drop")
      
      
    } else if (method == "GEE") {
      eta <- data1 %>% mutate(eta1= predict(model, newdata = mutate(data1, A = 1), type = "response")) %>%
        mutate(eta0 = predict(model, newdata = mutate(data1, A = 0), type = "response")) %>%
        dplyr::select(cluster, eta1, eta0) %>%
        group_by(cluster) %>%
        summarise_all(mean) %>%
        as.data.frame()
    } else if (method == "GLMM") {
      #pred1 <- as.matrix(cbind(rep(1,nrow(data1)),rep(1,nrow(data1)),data1[,c(grep("^X\\d+c$", names(data1)), grep("^X\\d+b$", names(data1)), grep("^H\\d+$", names(data1)))]))%*%as.vector(fixef(model))
      #pred0 <- as.matrix(cbind(rep(1,nrow(data1)),rep(0,nrow(data1)),data1[,c(grep("^X\\d+c$", names(data1)), grep("^X\\d+b$", names(data1)), grep("^H\\d+$", names(data1)))]))%*%as.vector(fixef(model))
      pred1 <- predict(model, newdata = transform(data1, A = 1), re.form = NA, type = "link")
      pred0 <- predict(model, newdata = transform(data1, A = 0), re.form = NA, type = "link")
      
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
