simulate_data <- function(
    n               = 5000,
    treat_override  = c("random","all_treated","all_untreated"),
    ## --- NEW options ---------------------------------------------------
    np_hazard       = FALSE,      # non-proportional treatment effect?
    HR_early        = 1.50,       # HR (treated vs ctrl) 0-90 d if np_hazard
    HR_late         = 0.80,       # HR  > 90 d           if np_hazard
    dep_censor      = FALSE,      # informative admin censoring?
    censor_base     = 1/100,      # mean 100 d when dep_censor = FALSE
    seed            = 42)
{
  set.seed(seed)
  treat_override <- match.arg(treat_override)

  #### 1) Baseline covariates  (unchanged) ##############################
  age <- pmax(rnorm(n, 48, 13), 18)
  sex <- rbinom(n, 1, 0.58)
  race_levels <- c("white","black","hispanic","asian","other")
  race_probs  <- c(.48,.14,.06,.02,.30)
  race_raw    <- sample(race_levels, n, TRUE, race_probs)
  race_raw[sample.int(n, floor(.20*n))] <- NA           # 20 % missing
  region_levels <- c("NE","MW","S","W")
  region_probs  <- c(.20,.18,.37,.25)
  region_raw    <- sample(region_levels, n, TRUE, region_probs)

  # correlated covariates (as in your script) --------------------------
  logit_ment <- -1 + 0.01*(60 - age) + 0.2*(1-sex)
  ment_ill   <- rbinom(n,1, plogis(logit_ment))
  substance_abuse <- rbinom(n,1, plogis(-2 + 1*ment_ill))
  smoking     <- rbinom(n,1, plogis(-1.5 + .7*substance_abuse + .3*sex))
  bmi         <- rnorm(n, 28, 5)
  over_obesity<- as.integer(bmi > 30)
  diabetes    <- rbinom(n,1, plogis(-2  + .02*age + .5*over_obesity))
  hypertension<- rbinom(n,1, plogis(-1.5+ .03*age + .4*over_obesity))
  cirrhosis   <- rbinom(n,1, plogis(-2  + 1.2*substance_abuse))
  ckd         <- rbinom(n,1, plogis(-3  + .04*age + .5*hypertension))
  hiv         <- rbinom(n,1, .04)
  copd        <- rbinom(n,1, .20)
  cancer      <- rbinom(n,1, .07)
  baseline_gfr<- rnorm(n, 90, 15)

  #### 2) Treatment assignment ##########################################
  logit_treat <- -1.2 + .015*age + .2*sex + .4*cirrhosis + .5*ckd +
    .3*diabetes + .2*hypertension + .3*ment_ill +
    .2*substance_abuse + .4*(ckd*cirrhosis) + .01*baseline_gfr
  logit_treat[region_raw=="S"] <- logit_treat[region_raw=="S"] - .2
  logit_treat[region_raw=="W"] <- logit_treat[region_raw=="W"] + .1
  logit_treat[!is.na(race_raw) & race_raw=="black"]  <- logit_treat[!is.na(race_raw) & race_raw=="black"]  + .1
  logit_treat[!is.na(race_raw) & race_raw=="asian"]  <- logit_treat[!is.na(race_raw) & race_raw=="asian"]  + .2
  logit_treat[!is.na(race_raw) & race_raw=="other"]  <- logit_treat[!is.na(race_raw) & race_raw=="other"]  - .1

  if (treat_override=="random") {
    treatment <- rbinom(n,1, plogis(logit_treat))
  } else if (treat_override=="all_treated") {
    treatment <- rep(1,n)
  } else {
    treatment <- rep(0,n)
  }

  #### 3) Event-follow_time generation #########################################
  # baseline log-hazard (prognostic score)
  lp_outcome <- -2.5 + .03*age + .6*ckd + .4*cirrhosis + .2*ment_ill +
    .3*diabetes + .1*substance_abuse
  base_rate  <- 0.02*exp(lp_outcome)          # λ0i

  if (!np_hazard) {
    # --- proportional hazards ---------------------------------------
    rate_i <- base_rate * exp(.06*treatment + .04*(treatment*ckd))
    follow_time_to_event <- rexp(n, rate = rate_i)
  } else {
    # --- non-proportional hazards (piece-wise at 90 d) --------------
    rate1 <- base_rate * ifelse(treatment==1, HR_early, 1)
    rate2 <- base_rate * ifelse(treatment==1, HR_late,  1)
    rpexp_piece <- function(n, r1, r2, tau=90) {
      u <- runif(n); p1 <- 1-exp(-r1*tau); t <- numeric(n)
      early <- u<=p1
      t[early] <- -log(1-u[early])/r1[early]
      t[!early]<- tau - log((1-u[!early])/(1-p1[!early]))/r2[!early]
      t
    }
    follow_time_to_event <- rpexp_piece(n, rate1, rate2, 90)
  }

  #### 4) Administrative censoring #####################################
  if (!dep_censor) {
    censor_follow_time <- rexp(n, rate = censor_base)
  } else {
    # hazard ↑ with prognostic score + treatment
    cens_rate <- censor_base * exp(.4*lp_outcome + .3*treatment)
    censor_follow_time <- rexp(n, rate = cens_rate)
  }

  observed_follow_time <- pmin(follow_time_to_event, censor_follow_time)
  event         <- as.integer(follow_time_to_event <= censor_follow_time)

  #### 5) Data frame ####################################################
  df <- data.frame(age, sex, race=race_raw, region=region_raw,
                   ment_ill, substance_abuse, smoking,
                   over_obesity, diabetes, hypertension,
                   cirrhosis, ckd, hiv, copd, cancer,
                   baseline_gfr, bmi,
                   treatment, follow_time=observed_follow_time, event)

  #### 6) Optional missingness ##########################################
  if (TRUE) {
    p_region_miss <- ifelse(df$region=="MW", .01, .05)
    df$region[runif(n) < p_region_miss] <- NA
    df$ckd[ runif(n) < plogis(-3 + 1*df$ckd + .5*df$cirrhosis) ] <- NA
  }
  df
}


## proportional hazards, independent censoring  (original behaviour)
d0 <- simulate_data(n = 5000)


library(survival)
fit_full <- coxph(Surv(follow_time, event) ~ treatment, data = d0)
hr_full  <- exp(coef(fit_full))
hr_full


## non-proportional hazards *and* dependent censoring
d1 <- simulate_data(n = 5000, np_hazard = TRUE, dep_censor = TRUE)
write.csv(d1, file="data/sim_hcv_aki_nph.csv", row.names = FALSE)


fit_full <- coxph(Surv(follow_time, event) ~ treatment, data = d1)
hr_full  <- exp(coef(fit_full))
hr_full


## force everyone treated, keep NP hazard – for “truth” calculation
d_truth1 <- simulate_data(n = 1e5, treat_override = "all_treated",np_hazard = TRUE, dep_censor = F)
d_truth0 <- simulate_data(n = 1e5, treat_override = "all_untreated",np_hazard = TRUE, dep_censor = F)
mean(d_truth1$event)/mean(d_truth0$event)


