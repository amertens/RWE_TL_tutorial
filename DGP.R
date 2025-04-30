
library(tidyverse)


##############################################################
##  generate_hcv_data()  —  HCV-AKI plasmode generator      ##
##############################################################
#  Args
#   N              cohort size *after* baseline exclusions
#   max_follow     administrative window (days)
#   risk_window    grace period after switch (days)
#   h0             baseline hazard – set lower to reduce event prevalence
#   HR_early       hazard multiplier (SOF vs non)   0–90 d
#   HR_late        multiplier beyond 90 d
#   p_sof          target marginal prevalence of SOF (if treat_override = "simulate")
#   impute         logical; use missForest to impute region & CKD
#   treat_override     "simulate", "all_treated", or "all_control"
#   seed           RNG seed
#
#  Returns         tibble ready for analysis  (one row per subject)
generate_hcv_data <- function(
    ## core knobs --------------------------------------------------------
    N            = 125000,
    p_sof        = 0.36,
    h0           = 3.5e-4,
    HR_early     = 1.50,     # causal HR up to τ
    HR_late      = 0.70,     # causal HR after τ  (ignored if np_hazard = FALSE)
    tau          = 90,       # change-point (days)
    max_follow   = 180,
    risk_window  = 30,
    ## NEW options -------------------------------------------------------
    np_hazard    = FALSE,    # time-varying treatment effect?
    dep_censor   = FALSE,    # censoring depends on risk + treatment?
    censor_base  = 1/100,    # baseline censoring rate
    ## misc --------------------------------------------------------------
    treat_override = c("simulate","all_treated","all_control"),
    add_missing  = FALSE,
    impute       = FALSE,
    seed         = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  treat_override <- match.arg(treat_override)
  if (impute) requireNamespace("missForest")
  requireNamespace("tidyverse")

  ## 1  Demography -------------------------------------------------------
  raw <- tibble::tibble(
    id       = seq_len(N),
    age      = pmax(rnorm(N, 48, 13), 18),
    sex_male = rbinom(N, 1, 0.58),
    race     = sample(c("white","black","hispanic","asian","other"), N, TRUE,
                      prob = c(.48,.14,.06,.02,.30)),
    region   = sample(c("NE","MW","S","W"), N, TRUE,
                      prob = c(.20,.18,.37,.25)),
    enroll_days = rpois(N, 420)
  )

  ## 2  Clinical history & meds (unchanged) ------------------------------
  add_bin <- function(p) rbinom(N, 1, p)
  raw <- raw %>%
    mutate(
      ckd     = add_bin(.08),  prior_aki = add_bin(.05),
      heart_failure = add_bin(.07), sepsis = add_bin(.03),
      dehydration = add_bin(.06), obstruction = add_bin(.04),
      cirrhosis = add_bin(.18), portal_htn = add_bin(.04),
      esld   = add_bin(.02), hiv = add_bin(.04),
      diabetes = add_bin(.20), hypertension = add_bin(.45),
      bmi    = rnorm(N, 28, 5), overweight_obese = add_bin(.20),
      smoking = add_bin(.40), alcohol = add_bin(.18),
      substance_abuse = add_bin(.25), cancer = add_bin(.08),
      chemo  = add_bin(.01),
      nsaid  = add_bin(.25), acearb = add_bin(.30), diuretic = add_bin(.22),
      aminoglycoside = add_bin(.05), contrast = add_bin(.08),
      statin = add_bin(.15), aspirin = add_bin(.10),
      beta_blocker = add_bin(.14), ccb = add_bin(.16), art = add_bin(.05),
      prior_sof = add_bin(.05), prior_nonsof = add_bin(.05)
    )

  ## 3  Baseline exclusions --------------------------------------------
  cohort <- raw %>%
    filter(enroll_days >= 365,
           age >= 18,
           prior_aki == 0,
           !(prior_sof == 1 | prior_nonsof == 1))

  ## 4  Treatment assignment -------------------------------------------
  if (treat_override == "simulate") {
    lp <- with(cohort,
               0.015*age +
                 0.30*cirrhosis + 0.35*portal_htn + 0.25*ckd +
                 0.20*hiv + 0.15*substance_abuse + 0.10*diabetes +
                 0.05*hypertension -0.10*cancer -0.05*overweight_obese +
                 if_else(region=="W", 0.05, if_else(region=="S", -0.05, 0)) +
                 case_when(race=="black" ~ 0.05,
                           race=="asian" ~ 0.10,
                           race=="other" ~ -0.05,
                           TRUE ~ 0) +
                 rnorm(nrow(cohort), 0, 0.6))
    alpha0 <- qlogis(p_sof) - mean(lp)
    p_trt  <- plogis(alpha0 + lp) |> pmin(0.95) |> pmax(0.05)
    cohort$treatment <- rbinom(nrow(cohort), 1, p_trt)
  } else {
    cohort$treatment <- ifelse(treat_override=="all_treated",1L,0L)
  }

  ## 5  Individual baseline hazard -------------------------------------
  lp_out <- with(cohort,
                 -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis +
                   0.3*heart_failure + 0.25*nsaid + 0.20*contrast)
  base_rate <- h0 * exp(lp_out)

  ## 6  Event times -----------------------------------------------------
  if (!np_hazard) {
    # proportional hazards
    rate  <- base_rate * ifelse(cohort$treatment==1, HR_early, 1)
    cohort$event_time <- rexp(nrow(cohort), rate = rate)
  } else {
    # piece-wise HR (early vs late)
    rpexp_piece <- function(n, r1, r2, tau){
      u  <- runif(n); p1 <- 1-exp(-r1*tau); t <- numeric(n)
      early <- u <= p1
      t[early]  <- -log(1-u[early]) / r1[early]
      t[!early] <- tau - log((1-u[!early])/(1-p1[!early])) / r2[!early]
      t
    }
    r1 <- base_rate * ifelse(cohort$treatment==1, HR_early, 1)
    r2 <- base_rate * ifelse(cohort$treatment==1, HR_late,  1)
    cohort$event_time <- rpexp_piece(nrow(cohort), r1, r2, tau)
  }

  ## 7  Administrative censoring (independent or informative) ----------
  if (!dep_censor) {
    censor_admin <- rexp(nrow(cohort), rate = censor_base)
  } else {
    cens_rate <- censor_base * exp(0.4*lp_out + 0.3*cohort$treatment)
    censor_admin <- rexp(nrow(cohort), rate = cens_rate)
  }
  cohort$censor_admin <- pmin(censor_admin, max_follow)

  ## switch censoring stays as before
  cohort$tx_days <- ifelse(cohort$treatment==1, rpois(nrow(cohort),84),
                           rpois(nrow(cohort),70))
  cohort$switch  <- rbinom(nrow(cohort), 1, .03)
  cohort$censor_switch <- ifelse(cohort$switch==1,
                                 cohort$tx_days + risk_window,
                                 max_follow)

  cohort$follow_time <- pmin(cohort$event_time,
                             cohort$censor_admin,
                             cohort$censor_switch)
  cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

  ## 8  Analysis data ---------------------------------------------------
  ana <- cohort %>%
    dplyr::select(-enroll_days, -prior_aki, -prior_sof, -prior_nonsof,
                  -tx_days, -event_time, -censor_admin, -censor_switch)

  ## 9  Optional missingness & imputation -------------------------------
  if (add_missing) {
    ana$region[runif(nrow(ana)) < .05] <- NA
    ana$ckd[   runif(nrow(ana)) < .10] <- NA
    if (impute) {
      imp_vars <- c("age","race","region","ckd","cirrhosis","hiv",
                    "diabetes","hypertension","bmi")
      imp_in <- ana[,imp_vars] %>% mutate(across(c(race,region), as.factor))
      ana[,imp_vars] <- missForest::missForest(as.data.frame(imp_in),
                                               verbose = FALSE)$ximp
    }
  }
  ana
}
