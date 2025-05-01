library(tidyverse)

N             = 125000
p_sof         = 0.36
h0            = 3.5e-4
HR_early      = 1.50
HR_late       = 0.70
tau           = 90
max_follow    = 180
risk_window   = 30
## NEW options -------------------------------------------------------
np_hazard     = FALSE
dep_censor    = FALSE
censor_base   = 1/100
complexity    = FALSE
## misc --------------------------------------------------------------
treat_override= "simulate"
add_missing   = FALSE
impute        = FALSE
seed          = NULL

generate_hcv_data <- function(
    ## core knobs --------------------------------------------------------
    N             = 125000,
    p_sof         = 0.36,
    h0            = 3.5e-4,
    HR_early      = 1.50,
    HR_late       = 0.70,
    tau           = 90,
    max_follow    = 180,
    risk_window   = 30,
    ## NEW options -------------------------------------------------------
    np_hazard     = FALSE,
    dep_censor    = FALSE,
    censor_base   = 1/100,
    complexity    = FALSE,
    ## misc --------------------------------------------------------------
    treat_override= c("simulate","all_treated","all_control"),
    add_missing   = FALSE,
    impute        = FALSE,
    seed          = NULL
){
  if (!is.null(seed)) set.seed(seed)
  treat_override <- match.arg(treat_override)
  if (impute) requireNamespace("missForest")
  requireNamespace("tidyverse")

  ## 1  Demography -------------------------------------------------------
  raw <- tibble(
    id           = seq_len(N),
    age          = pmax(rnorm(N, 48, 13), 18),
    sex_male     = rbinom(N, 1, 0.58),
    race         = sample(c("white","black","hispanic","asian","other"), N, TRUE,
                          prob = c(.48,.14,.06,.02,.30)),
    region       = sample(c("NE","MW","S","W"), N, TRUE,
                          prob = c(.20,.18,.37,.25)),
    enroll_days  = rpois(N, 420),
    Nobs=N
  )

  ## 2  Clinical history & meds -----------------------------------------
  add_bin <- function(p) rbinom(N, 1, p)
  raw <- raw %>%
    mutate(
      ckd     = add_bin(.08),  prior_aki = add_bin(.05),
      heart_failure = add_bin(.07), sepsis = add_bin(.03),
      dehydration   = add_bin(.06), obstruction = add_bin(.04),
      cirrhosis      = add_bin(.18), portal_htn = add_bin(.04),
      esld   = add_bin(.02), hiv   = add_bin(.04),
      diabetes = add_bin(.20), hypertension = add_bin(.45),
      bmi       = rnorm(N, 28, 5), overweight_obese = add_bin(.20),
      smoking   = add_bin(.40), alcohol = add_bin(.18),
      substance_abuse = add_bin(.25), cancer = add_bin(.08),
      chemo      = add_bin(.01),
      nsaid      = add_bin(.25), acearb = add_bin(.30), diuretic = add_bin(.22),
      aminoglycoside = add_bin(.05), contrast    = add_bin(.08),
      statin     = add_bin(.15), aspirin = add_bin(.10),
      beta_blocker = add_bin(.14), ccb    = add_bin(.16), art = add_bin(.05),
      prior_sof    = add_bin(.05), prior_nonsof = add_bin(.05)
    )

  ## 3  Baseline exclusions --------------------------------------------
  cohort <- raw %>%
    filter(enroll_days >= 365, age >= 18, prior_aki == 0,
           !(prior_sof == 1 | prior_nonsof == 1))

  ## 4  Treatment assignment -------------------------------------------
  if (treat_override=="simulate") {
    ## linear lp for simple case
    lp0 <- with(cohort,
                0.015*age + 0.30*cirrhosis + 0.25*ckd + 0.15*hiv + 0.10*diabetes -
                  0.10*cancer + rnorm(Nobs,0,0.6)
    )
    if (complexity) {
      ## add nonlinear & interaction terms
      lp0 <- lp0 +
        0.02*(cohort$bmi^2)/100 -
        0.3*sin(0.1*cohort$bmi) +
        0.5*(cohort$age/50)^3 +
        1.5*cohort$ckd*cohort$cancer +    # interaction
        0.8*cohort$hiv * log1p(cohort$age)
    }
    alpha0 <- qlogis(p_sof) - mean(lp0)
    p_trt  <- plogis(alpha0 + lp0) |> pmin(0.95) |> pmax(0.05)
    cohort$treatment <- rbinom(nrow(cohort),1,p_trt)
  } else {
    cohort$treatment <- ifelse(treat_override=="all_treated",1L,0L)
  }

  ## 5  Individual baseline hazard -------------------------------------
  if (!complexity) {
    lp_out <- with(cohort,
                   -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis +
                     0.3*heart_failure + 0.25*nsaid + 0.20*contrast
    )
  } else {
    ## add nonlinear hazard terms + interactions
    lp_out <- with(cohort,
                   -2.8 +
                     0.03*age + 0.0005*age^2 +
                     0.7*ckd + 0.5*cirrhosis +
                     0.02*(bmi^2)/100 -
                     0.3*sin(0.1*bmi) +
                     0.4*heart_failure*acearb +     # interaction
                     0.6*nsaid*cohort$treatment +
                     0.3*contrast*log1p(age)
    )
  }
  base_rate <- h0 * exp(lp_out)

  ## 6  Event times -----------------------------------------------------
  if (!np_hazard) {
    rate <- base_rate * ifelse(cohort$treatment==1,HR_early,1)
    cohort$event_time <- rexp(nrow(cohort), rate=rate)
  } else {
    ## piecewise hazard as before...
    rpexp_piece <- function(n,r1,r2,tau){
      u <- runif(n); p1 <- 1-exp(-r1*tau); t <- numeric(n)
      e <- u <= p1
      t[e]  <- -log(1-u[e]) / r1[e]
      t[!e] <- tau - log((1-u[!e])/(1-p1[!e]))/r2[!e]
      t
    }
    r1 <- base_rate * ifelse(cohort$treatment==1,HR_early,1)
    r2 <- base_rate * ifelse(cohort$treatment==1,HR_late,1)
    cohort$event_time <- rpexp_piece(nrow(cohort),r1,r2,tau)
  }

  ## 7  Admin censoring -------------------------------------------------
  if (!dep_censor) {
    censor_admin <- rexp(nrow(cohort), rate=censor_base)
  } else {
    c_rate <- censor_base * exp(0.4*lp_out + 0.3*cohort$treatment)
    censor_admin <- rexp(nrow(cohort), rate=c_rate)
  }

  cohort$censor_admin <- pmin(censor_admin, max_follow)

  cohort$tx_days  <- ifelse(cohort$treatment==1,
                            rpois(N,84), rpois(N,70))
  cohort$switch   <- rbinom(nrow(cohort),1,0.03)
  cohort$censor_switch <- ifelse(cohort$switch==1,
                                 cohort$tx_days + risk_window, max_follow)

  cohort$follow_time <- pmin(cohort$event_time,
                             cohort$censor_admin,
                             cohort$censor_switch)
  cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

  ## 8  Analysis dataset ------------------------------------------------
  ana <- cohort %>%
    select(-enroll_days,-prior_aki,-prior_sof,-prior_nonsof,
           -tx_days,-event_time,-censor_admin,-censor_switch)

  ## 9  Missingness & imputation ----------------------------------------
  if (add_missing) {
    ana$region[sample(Nobs, 0.05*Nobs)] <- NA
    ana$ckd[sample(Nobs, 0.10*Nobs)] <- NA
    if (impute) {
      imp_vars <- c("age","race","region","ckd","cirrhosis","hiv",
                    "diabetes","hypertension","bmi")
      imp_in   <- ana %>% select(all_of(imp_vars)) %>%
        mutate(across(c(race,region), as.factor))
      ana[,imp_vars] <- missForest::missForest(as.data.frame(imp_in),
                                               verbose=FALSE)$ximp
    }
  }

  return(ana)
}


df <- generate_hcv_data(    np_hazard     = FALSE,
                            dep_censor    = FALSE,
                            complexity    = FALSE,
                            add_missing   = FALSE)

write.csv(df,here::here("data/sim_hcv_aki.csv"))
