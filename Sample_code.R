## This file contains functions calculating point estimates, associated 95% CIs,
## and p-values for testing the non-informative cluster sizes for the methods 
## described in the manuscript titled ‘Model-Robust Standardization in Cluster-Randomized Trials.’

#Function calculating point estimates
robust_CRT_point <- function(m_formula, df, clus_id, trt, prob, out_model, scale){

  #Input:
  #m_formula: formula of outcome mean model
  #df: data frame where the categorical variables should be already converted to dummy variables
  #clus_id: a character of the variable name of the cluster id
  #trt: a character of the variable name of the treatment assignment per cluster
  #prob: a vector of probabilities of treatment assignment per cluster conditional on covariates
  #out_model: specifications of outcome mean models; 
  #potential values are: (i) continuous outcomes: 'LM' (linear model on cluster-level means), 'LMM' (linear mixed model on individual-level data), 'GEE_EX_CON' (generalized estimating equations with exchangeable correlation), 'GEE_IND_CON' (generalized estimating equations with independent correlation);
  #(ii) and binary outcomes: 'GLMM' (generalized linear mixed model on individual-level data), 'GEE_EX_BIN' (generalized estimating equations with exchangeable correlation), 'GEE_IND_BIN' (generalized estimating equations with independent correlation)
  #scale: risk differences ('RD'), relative risks ('RR') and odds ratios ('OR')
  
  require(dplyr)
  require(geepack) 
  require(lme4) 
  require(nlme)
  
  # Validate inputs
  stopifnot(is.data.frame(df), is.character(clus_id), is.character(outcome), is.character(trt))
  stopifnot(is.vector(prob), is.character(out_model))
  
  # Extract all variable names from the formula
  m_formula <- Y~A+X1+X2+H1+H2
  all_vars <- all.vars(formula(m_formula))
  #outcome
  outcome <- all_vars[1]
  #covariates
  cov <- setdiff(all_vars, c(outcome, trt))

  # Check if a variable is cluster-level
  is_cluster_level <- function(var, group) {
    tapply_result <- tapply(var, group, function(x) length(unique(x)))
    all(tapply_result == 1)  # Check if all groups have unique length 1
  }
  clus_cov <- cov[sapply(cov, function(covariate) is_cluster_level(df[[covariate]], df[[clus_id]]))]
  ind_cov <- setdiff(cov, clus_cov)
  
  # Create new data set for analysis
  df1 <- df %>% select(all_of(c(clus_id, trt, outcome, clus_cov, ind_cov))) %>%
    rename(clus_id = !!clus_id, A = !!trt, Y = !!outcome) %>%
    arrange(clus_id)
  df1$prob <- prob
  ind_cov_names <- paste0("X", seq_along(ind_cov))
  clus_cov_names <- paste0("H", seq_along(clus_cov))
  names(df1)[4:(3 + length(ind_cov))] <- ind_cov_names
  names(df1)[(4 + length(ind_cov)):(3 + length(ind_cov) + length(clus_cov))] <- clus_cov_names
  
  # Calculate cluster sizes
  cluster_sizes <- df1 %>% group_by(clus_id) %>% summarise(N = n())
  df1 <- left_join(df1, cluster_sizes, by = "clus_id")
  
  # Genrate a data frame of cluster means 
  df_clus <- as.data.frame(df1 %>% group_by(clus_id) %>% summarise(across(everything(), mean)))
  for (col in ind_cov_names) {
    df1 <- df1 %>%
      mutate(
        !!paste0(col, "b") := rep(df_clus[[col]], df_clus$N),
        !!paste0(col, "c") := .data[[col]] - .data[[paste0(col, "b")]]
      )
  }
  
  # Build formula for individual-level models
  formulai <- paste("Y ~ A +",
    paste(paste0(ind_cov_names, "c"), collapse = " + "), " + ",
    paste(paste0(ind_cov_names, "b"), collapse = " + "), " + ",
    paste(clus_cov_names, collapse = " + "), " + N"
  )
  # Build formula for cluster-level models
  formulac <- paste("Y ~ A +",
    paste(ind_cov_names, collapse = " + "), " + ",
    paste(clus_cov_names, collapse = " + "), " + N"
  )
  
  # Fit models
  model <- switch(out_model,
                  "LM" = try(lm(formulac, data = df_clus), silent = T),
                  "LMM" = try(lme(as.formula(formulai), random = ~ 1 | clus_id, data = df1), silent = T),
                  "GEE_EX_CON" = try(geeglm(as.formula(formulai), id = clus_id, data = df1, corstr = "exchangeable"), silent=T),
                  "GEE_IND_CON" = try(geeglm(as.formula(formulai), id = clus_id, data = df1, corstr = "independence"), silent=T),
                  "GLMM" =  try(glmer(paste(formulai, "+ (1 | clus_id)"), data = df1, binomial(link = "logit"), nAGQ = 0)),
                  "GEE_EX_BIN" = try(geeglm(as.formula(formulai), id = clus_id, data = df1, corstr = "exchangeable", family = binomial("logit"))),
                  "GEE_IND_BIN" = try(glm(as.formula(formulai), data = df1, family = binomial("logit"))),
                  stop("Invalid model type specified.")
  )
  
  if (out_model == "LM"){
    #eta: vector containing the predicted outcome in two arms 
    eta <- df_clus %>% mutate(eta1= predict(model, newdata = mutate(df_clus, A = 1), type = "response")) %>%
      mutate(eta0 = predict(model, newdata = mutate(df_clus, A = 0), type = "response")) %>%
      dplyr::select(clus_id, eta1, eta0) %>%
      group_by(clus_id) %>%
      as.data.frame()
  } else if (out_model == "LMM") {
    eta <- df_clus %>% mutate(eta1= as.matrix(cbind(rep(1,nrow(df_clus)),rep(1,nrow(df_clus)),df_clus[,c(grep("^X\\d$", names(df_clus)), grep("^H", names(df_clus)), grep("^N$", names(df_clus)))]))%*%as.vector(fixef(model)[-c(grep("^X\\d+c$", names(fixef(model))))])) %>%
      mutate(eta0 = as.matrix(cbind(rep(1,nrow(df_clus)),rep(0,nrow(df_clus)),df_clus[,c(grep("^X\\d$", names(df_clus)), grep("^H", names(df_clus)), grep("^N$", names(df_clus)))]))%*%as.vector(fixef(model)[-c(grep("^X\\d+c$", names(fixef(model))))])) %>%
      as.data.frame()
  } else if (out_model == "GEE_EX_CON" | out_model == "GEE_IND_CON" | out_model == "GEE_EX_BIN" | out_model == "GEE_IND_BIN") {
    eta <- df1 %>% mutate(eta1= predict(model, newdata = mutate(df1, A = 1), type = "response")) %>%
      mutate(eta0 = predict(model, newdata = mutate(df1, A = 0), type = "response")) %>%
      dplyr::select(clus_id, eta1, eta0) %>%
      group_by(clus_id) %>%
      summarise_all(mean) %>%
      as.data.frame()
  } else if (out_model == "GLMM") {
    #pred1 <- as.matrix(cbind(rep(1,nrow(df1)),rep(1,nrow(df1)),df1[,c(grep("^X\\d$", names(df1)), grep("^H\\d+$", names(df1)), grep("^N$", names(df1)))]))%*%as.vector(fixef(model))
    #pred0 <- as.matrix(cbind(rep(1,nrow(df1)),rep(0,nrow(df1)),df1[,c(grep("^X\\d$", names(df1)), grep("^H\\d+$", names(df1)), grep("^N$", names(df1)))]))%*%as.vector(fixef(model))
    pred1 <- as.matrix(cbind(rep(1,nrow(df1)),rep(1,nrow(df1)),df1[,c(grep("^X\\d+c$", names(df1)), grep("^X\\d+b$", names(df1)), grep("^H\\d+$", names(df1)), grep("^N$", names(df1)))]))%*%as.vector(fixef(model))
    pred0 <- as.matrix(cbind(rep(1,nrow(df1)),rep(0,nrow(df1)),df1[,c(grep("^X\\d+c$", names(df1)), grep("^X\\d+b$", names(df1)), grep("^H\\d+$", names(df1)), grep("^N$", names(df1)))]))%*%as.vector(fixef(model))
    #pi: mathematical constant
    Pi <- 3.141592653589793
    eta1 <- df1 %>% mutate(eta1 = (exp(pred1/sqrt(3*as.numeric(VarCorr(model))/Pi^2+1))/(1+exp(pred1/sqrt(3*as.numeric(VarCorr(model))/Pi^2+1))))) %>%
      group_by(clus_id) %>%
      summarise_all(mean) %>%
      as.data.frame() 
    eta0 <-  df1 %>% mutate(eta0 = (exp(pred0/sqrt(3*as.numeric(VarCorr(model))/Pi^2+1))/(1+exp(pred0/sqrt(3*as.numeric(VarCorr(model))/Pi^2+1))))) %>%
      group_by(clus_id) %>%
      summarise_all(mean) %>%
      as.data.frame() 
    eta <- df_clus %>% mutate(eta1 = eta1[,ncol(eta1)]) %>%
      mutate(eta0 = eta0[,ncol(eta0)]) %>%
      dplyr::select(clus_id, eta1, eta0) %>%
      as.data.frame()
  } 
  
  #####################Point estimates using our proposed methods###############
  mu_C1 <- df_clus$A / df_clus$prob * (df_clus$Y - eta$eta1) + eta$eta1
  mu_C0 <- (1-df_clus$A) / (1-df_clus$prob) * (df_clus$Y - eta$eta0) + eta$eta0
  cate <- switch(scale,
                 "RD" = mean(mu_C1) - mean(mu_C0),
                 "RR" = mean(mu_C1)/mean(mu_C0),
                 "OR" =  mean(mu_C1)/(1-mean(mu_C1))/mean(mu_C0)*(1-mean(mu_C0)),
                 stop("Invalid scales specified."))
  mu_I1 <- mean(df_clus$N * mu_C1)/mean(df_clus$N)
  mu_I0 <- mean(df_clus$N * mu_C0)/mean(df_clus$N)
  iate <- switch(scale,
                 "RD" = mean(mu_I1) - mean(mu_I0),
                 "RR" = mean(mu_I1)/mean(mu_I0),
                 "OR" =  mean(mu_I1)/(1-mean(mu_I1))/mean(mu_I0)*(1-mean(mu_I0)),
                 stop("Invalid scales specified."))
  
  test_NICS <-  cate - iate
  return(list(df1,
              df_clus,
              c(cate,iate,test_NICS)))
}

robust_CRT <- function(m_formula, df, clus_id, trt, prob, out_model, scale, jack = 1, alpha){
  
  #Input:
  #m_formula: formula of outcome mean model
  #df: data frame where the categorical variables should be already converted to dummy variables
  #clus_id: name of the variable representing cluster id
  #trt: name of the variable representing treatment assignment per cluster
  #prob: a vector of probabilities of conditional treatment assignment per cluster
  #out_model: potential choices of outcome regression; 
  #potential values include continuous outcomes: 'LM' (linear model on cluster-level means), 'LMM' (linear mixed model on individual-level data), 'GEE_EX_CON' (generalized estimating equations with exchangeable correlation), 'GEE_IND_CON' (generalized estimating equations with independent correlation),
  #and binary outcomes: 'GLMM' (generalized linear mixed model on individual-level data), 'GEE_EX_BIN' (generalized estimating equations with exchangeable correlation), 'GEE_IND_BIN' (generalized estimating equations with independent correlation)
  #scale: risk differences ('RD'), relative risks ('RR') and odds ratios ('OR')
  #jack: types of jackknife standard error estimates, potential values are 1, 2 and 3. Type 1 is the standard jackknife and type 3 is recommended for small number of clusters
  #alpha: type-I error rate
  
  temp <- robust_CRT_point(m_formula, df, clus_id, trt, prob, out_model, scale)
  df1 <- temp[[1]]
  df_clus <- temp[[2]]
  m <- nrow(df_clus)
  pes <- temp[[3]]
  cluster_names <- unique(df_clus$clus_id)
  point_est_jack <- matrix(NA, nrow = m, ncol = 3)
  for (i in cluster_names){
    df_jack <- df %>% filter(df[[clus_id]] != i)
    point_est_jack[which(cluster_names==i),] <- robust_CRT_point(m_formula, df_jack, clus_id, trt, prob[which(df[[clus_id]]!=i)], out_model, scale)[[3]]
  }
  # Jackknife standard error estimation
  if (jack == 1) {
    # Implement standard jackknife
    jackse <- sqrt(apply(point_est_jack, 2, var)*(m-1)^2/m)
  } else if (jack == 2) {
    # Implement alternative jackknife method 2
    rep_point <- as.matrix(do.call(rbind, replicate(m, pes, simplify = FALSE)))
    jackse <-  as.numeric(sqrt(colSums((point_est_jack-rep_point)^2)*(m-1)/m))
  } else if (jack == 3) {
    # Implement alternative jackknife method 3 
    rep_point <- as.matrix(do.call(rbind, replicate(m, pes, simplify = FALSE)))
    jackse2 <-  as.numeric(sqrt(colSums((point_est_jack-rep_point)^2)*(m-1)/m))
    jackse <- as.numeric(jackse2*sqrt(m/(m-1)))
  } else {
    stop("Invalid jackknife type specified.")
  }
  thres <- qt(p = 1-alpha/2, df = m-1)
  ciu <- pes+thres*jackse
  cil <- pes-thres*jackse
  table <- cbind(pes, jackse, cil, ciu)[-3,]
  rownames(table) <-  c("cATE",
                        "iATE")
  colnames(table) <- c("Estimate",
                       "Std, Error",
                       "CI lower",
                       "CI upper")
  table <- as.data.frame(round(table, 3))
  #test statistic for NICS
  test_sta <- pes[3]/jackse[3]
  p_val <- min((1-pt(test_sta, df = m-1, ncp = 0)),pt(test_sta, df = m-1, ncp = 0))*2
  table$`z value` <- c(round(test_sta, 3), NA) 
  table$`p value` <- c(round(p_val, 3), NA)
  table[is.na(table)] <- ""
  
  return(table)
}
out_model <- "LM"
robust_CRT(m_formula, df, clus_id, trt, prob, out_model, scale, jack = 1, 0.05)

