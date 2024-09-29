source("robust_CRT.R")
source("CRT_point.R")

library(rio)


require(dplyr)
require(geepack) 
require(lme4) 
require(nlme)



# Load the RData file and rename the object as 'cont_dt'
cont_dt <- import("cont_data_example.RData")
bin_dt <- import("bin_data_example.RData")

# Prepare data and calculate probabilities
df2 <- cont_dt %>%
  group_by(cluster_id) %>%
  mutate(first_A = first(A)) %>%
  ungroup() %>%
  mutate(prob_A_1 = mean(first_A == 1, na.rm = TRUE),  # Proportion A=1
         prob_A_0 = mean(first_A == 0, na.rm = TRUE)) %>%
  mutate(assigned_value = ifelse(A == 1, prob_A_1, prob_A_0))

prob <- df2$assigned_value

# Example usage of robust_CRT
example1 <- robust_CRT(
  formula = Y ~ X1 + X2 + cluster(H1 + H2 + N),
  data = cont_dt,
  clus_id = "cluster_id",
  trt = "A",
  prob = prob,
  method = "LM",
  corstr = NA,
  scale = "RR"
)



# Prepare data and calculate probabilities
df2 <- bin_dt %>%
  group_by(cluster_id) %>%
  mutate(first_A = first(A)) %>%
  ungroup() %>%
  mutate(prob_A_1 = mean(first_A == 1, na.rm = TRUE),  # Proportion A=1
         prob_A_0 = mean(first_A == 0, na.rm = TRUE)) %>%
  mutate(assigned_value = ifelse(A == 1, prob_A_1, prob_A_0))

prob <- df2$assigned_value


example2 <- robust_CRT(
  formula = Y ~ X1 + X2 + cluster(H1 + H2 + N),
  data = bin_dt,
  clus_id = "cluster_id",
  trt = "A",
  prob = prob,
  family = binomial("logit"),
  method = "GLMM",
  corstr = "exchangeable",
  scale = "RD"
)

