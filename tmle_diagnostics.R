############################################################
##  TMLE WORKFLOW  –  SOF vs non-SOF and AKI (simulated)  ##
############################################################
## Author: <your name>    Date: <today>

## --------------------------- 0. Setup ---------------------------

## Install packages the first time you run the script
pkgs <- c("tidyverse", "data.table", "tlverse",  # tmle3, sl3, etc.
          "SuperLearner", "glmnet", "xgboost", "randomForest",
          "cobalt", "EValue", "ggplot2", "cowplot")
lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, dependencies = TRUE))

library(tidyverse)
library(data.table)
library(tlverse)        # loads sl3 + tmle3
library(SuperLearner)   # to build custom SL libraries if desired
library(cobalt)         # balance diagnostics
library(EValue)
library(tmle)

set.seed(2025)

## --------------------------- 1. Generate data with poor balance ---------------------------
source("DGP.R")  # loads generate_hcv_data()

## create a data set where treatment **ignores** confounding structure
df_bad <- generate_hcv_data(N = 3000,
                            treat_override = "simulate",
                            seed = 101)

## Now *scramble* treatment to induce extreme imbalance
df_bad$treatment <- rbinom(nrow(df_bad), 1, 0.80)    # 80 % get SOF regardless of X

## Quick look
table(df_bad$treatment)

## --------------------------- 2. Declare Nodes for tmle3 ---------------------------
## tmle3 uses a list that maps variables to "roles"
node_list <- list(
  W = setdiff(names(df_bad), c("treatment", "event", "follow_time")),
  A = "treatment",
  Y = "event",
  Ttilde = "follow_time"  # required for survival TMLE
)

## --------------------------- 3. Specify the Super Learner ---------------------------
## A small, diverse library – feel free to extend
sl_lib <- c(
  Lrnr_glm$new(),                         # main-terms logistic regression
  Lrnr_glmnet$new(alpha = 0),             # ridge
  Lrnr_glmnet$new(alpha = 1),             # lasso
  Lrnr_xgboost$new(nrounds = 200,
                   max_depth = 3,
                   eta = 0.05,
                   subsample = 0.8),
  Lrnr_randomForest$new(num.trees = 400)
)

## --------------------------- 4. Fit TMLE (survival, risk difference at 180 days) ---------
## tmle3_helper_survival() sets up the specification
t_max   <- 180
tmle_sp <- tmle3_Spec_survival(time = t_max,
                               contrast = list(A = 1, B = 0),
                               censoring_node = NULL,          # no extra censor node
                               outcome_type = "binary",
                               marginal = "RD")                # risk difference

## Build task, learner list, and fit
task  <- tmle_sp$make_tmle_task(df_bad, node_list)
learner_outcome <- make_learner(Stack, sl_lib)
learner_treat   <- make_learner(Stack, sl_lib)

learner_list <- list(Y = learner_outcome,
                     A = learner_treat)

fit <- tmle3(tmle_sp, task, learner_list)

## --------------------------- 5. Results ---------------------------
tmle_est <- fit$summary
print(tmle_est)

## --------------------------- 6. Diagnostics ------------------------------------------------

### 6a. PROPENSITY-SCORE OVERLAP
ps <- fit$likelihood$get_likelihoods(task, "A")[["A"]]  # predicted e(X)
ps_df <- data.frame(ps = ps,
                    trt = factor(df_bad$treatment))
g1 <- ggplot(ps_df, aes(ps, fill = trt)) +
  geom_histogram(bins = 40, position = "identity",
                 alpha = 0.4, colour = "black") +
  labs(title = "Propensity-score overlap", x = "e(X)")

### 6b. WEIGHTED STANDARDISED MEAN DIFFERENCES
## Use IPTW weights from the fitted propensity score
wt     <- ifelse(df_bad$treatment == 1, 1/ps, 1/(1-ps))
bal.tab <- bal.tab(x = df_bad[node_list$W],
                   treat = df_bad$treatment,
                   weights = wt,
                   method = "weighting",
                   estimand = "ATE")
print(bal.tab)
g2 <- love.plot(bal.tab, threshold = .1) +
  labs(title = "Weighted covariate balance (SMD)")

### 6c. CROSS-VALIDATED LOG-LOSS OF THE OUTCOME SL
cv_risk <- learner_outcome$cv_risk(task, loss_function = loss_loglik_binomial)
print(cv_risk)

### 6d. TMLE INFLUENCE-CURVE DISTRIBUTION
ic <- fit$IC
g3 <- ggplot(data.frame(ic = ic), aes(ic)) +
  geom_histogram(bins = 50, colour = "black", fill = "steelblue", alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "IC distribution (should be approx. mean 0)")

### 6e. E-VALUE
risk_rd <- tmle_est$tmle_est
se_rd   <- tmle_est$se
rr_est  <- 1 + risk_rd / mean(df_bad$event)  # crude transform
e_val   <- evalue(est = rr_est,
                  lo = rr_est - 1.96*se_rd,
                  hi = rr_est + 1.96*se_rd,
                  true = 1,
                  type = "RR")
print(e_val)

### 6f. Display plots side-by-side
cowplot::plot_grid(g1, g2, g3, nrow = 1)



# After tmle3 fit
Y_star <- fit$updates$Y_star  # targeted prediction
resid  <- df$event - Y_star

# Stack covariates
Xmat   <- as.matrix(df[ , node_list$W])

# Regress residuals on all covariates
lm_fit <- lm(resid ~ Xmat)

# Global F-test
anova(lm_fit)
