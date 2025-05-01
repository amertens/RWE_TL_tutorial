
library(tidyverse)
library(survivalSL)
library(survtmle)
library(here)
source(here("DGP.R"))

# 1) Simulate the data
set.seed(123)
dat <- generate_hcv_data(
  N             = 10000,         # sample size
  p_sof         = 0.36,          # marginal SOF prevalence
  h0            = 3.5e-4,        # baseline hazard
  HR_early      = 1.50,          # early HR
  HR_late       = 0.70,          # late HR (ignored if np_hazard=FALSE)
  tau           = 90,            # change point
  max_follow    = 180,           # administrative censoring window
  risk_window   = 30,            # grace period after switch
  np_hazard     = FALSE,         # proportional hazards
  dep_censor    = FALSE,         # independent censoring
  censor_base   = 1/100,
  treat_override= "simulate"     # simulate treatment assignment
)

# keep only the columns we need
df <- dat %>%
  select(id, treatment, age, sex_male, ckd, cirrhosis,
         hiv, diabetes, hypertension, nsaid, contrast,
         follow_time, event)

# 2) Fit survivalSL ----------------------------------------------------------

# define a small SL library
sl_lib <- c("SL.coxph","SL.rfsrc")

# We’ll split into a “training” and “validation” set just to demonstrate
# how to get out-of-sample survival curves.
set.seed(456)
train_idx <- sample(nrow(df), size = 0.7*nrow(df))
train_data <- df[train_idx,]
valid_data <- df[-train_idx,]

sl_fit <- survivalSL(
  methods=c("LIB_COXen", "LIB_AFTgamma", "LIB_PHexponential"),
  metric     = "ci",               # concordance index
  data       = as.data.frame(train_data),
  times      = "follow_time",      # time variable
  failures   = "event",            # 1=event, 0=censor
  cov.quanti = c("age"),           # continuous covariates
  cov.quali  = c("treatment","sex_male","ckd"), # binary/factor covariates
  progress   = TRUE
)

# Summarize predictive performance
summary(sl_fit, digits = 3)

# Extract the SuperLearner‐aggregated survival curves on the validation set
surv_curves <- predict(sl_fit, newdata = valid_data)
# surv_curves is a list containing predicted survival probabilities over time

# For example, plot the mean survival curve by treatment arm at a few timepoints:
time_grid <- seq(0,180,by=30)
pred_surv <- data.frame(
  time      = rep(time_grid, 2),
  surv_prob = c(
    colMeans(surv_curves$`1`[ , as.character(time_grid)]),  # SOF arm
    colMeans(surv_curves$`0`[ , as.character(time_grid)])   # non‐SOF arm
  ),
  treatment = rep(c("SOF","non‐SOF"), each = length(time_grid))
)
ggplot(pred_surv, aes(time, surv_prob, color=treatment)) +
  geom_line() +
  labs(x="Days since DAA start", y="Predicted survival",
       title="SuperLearner survival curves on validation set")

# 3) Fit survTMLE for an adjusted hazard ratio --------------------------------

# Prepare the adjustment matrix
W <- df %>%
  select(age, sex_male, ckd, cirrhosis, hiv, diabetes, hypertension,
         nsaid, contrast) %>%
  as.data.frame()

# Define Super Learner libraries for the nuisance functions
sl_lib_f <- c("SL.coxph","SL.rfsrc")  # failure model
sl_lib_c <- c("SL.coxph")             # censoring model

tmle_fit <- survtmle(
  ftime    = df$follow_time,     # follow-up or event/censoring time
  ftype    = df$event,           # 1=event(AKI), 0=censor
  trt      = df$treatment,       # 1=SOF, 0=non‐SOF
  adjustVars = W,                # baseline confounders
  t0        = 90,                # time point for HR estimation
  SL.ftime  = sl_lib_f,
  SL.ctime  = sl_lib_c,
  method    = "hazard"
)

# Inspect the HR estimate
print(tmle_fit$estimates)
# tmle_fit$estimates$psi.haz gives the TMLE hazard ratio at t0.

# 4) Compare with a naive Cox model -----------------------------------------

cox_fit <- coxph(
  Surv(follow_time, event) ~ treatment + age + sex_male + ckd +
    cirrhosis + hiv + diabetes + hypertension + nsaid + contrast,
  data = df
)
summary(cox_fit)

# The Cox coefficient on 'treatment' exp(coef) is the usual adjusted HR.

