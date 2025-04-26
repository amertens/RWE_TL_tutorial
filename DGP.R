simulate_data <- function(n = 5000,
                          treat_override = c("random","all_treated","all_untreated"),
                          # The rest of the parameters can remain hard-coded or made flexible.
                          seed = 42) {
  set.seed(seed)
  treat_override <- match.arg(treat_override)

  #### 1) Baseline Covariates ####
  age <- pmax(rnorm(n, mean=48, sd=13), 18)
  sex <- rbinom(n, 1, 0.58)

  race_levels <- c("white","black","hispanic","asian","other")
  race_probs  <- c(0.48, 0.14, 0.06, 0.02, 0.30)
  race_raw <- sample(race_levels, size=n, replace=TRUE, prob=race_probs)
  race_missing_ix <- sample.int(n, floor(0.20*n))
  race_raw[race_missing_ix] <- NA

  region_levels <- c("NE","MW","S","W")
  region_probs  <- c(0.20, 0.18, 0.37, 0.25)
  region_raw <- sample(region_levels, size=n, replace=TRUE, prob=region_probs)

  # Some correlated covariates (example)
  # mental illness depends on younger/female
  logit_ment <- -1 + 0.01*(60 - age) + 0.2*(1 - sex)
  prob_ment  <- plogis(logit_ment)
  ment_ill   <- rbinom(n,1,prob_ment)

  logit_subst <- -2 + 1.0*ment_ill
  substance_abuse <- rbinom(n,1, plogis(logit_subst))

  logit_smoke <- -1.5 + 0.7*substance_abuse + 0.3*sex
  smoking     <- rbinom(n,1, plogis(logit_smoke))

  # A few more
  bmi   <- rnorm(n, mean=28, sd=5)
  over_obesity <- as.integer(bmi > 30)

  logit_diab <- -2.0 + 0.02*age + 0.5*over_obesity
  diabetes   <- rbinom(n,1, plogis(logit_diab))

  logit_htn  <- -1.5 + 0.03*age + 0.4*over_obesity
  hypertension <- rbinom(n,1, plogis(logit_htn))

  logit_cirr <- -2 + 1.2*substance_abuse
  cirrhosis  <- rbinom(n,1, plogis(logit_cirr))

  logit_ckd  <- -3 + 0.04*age + 0.5*hypertension
  ckd        <- rbinom(n,1, plogis(logit_ckd))

  hiv <- rbinom(n,1,0.04)
  copd<- rbinom(n,1,0.20)
  cancer <- rbinom(n,1,0.07)
  baseline_gfr <- rnorm(n, mean=90, sd=15)

  #### 2) Treatment Assignment ####
  # If random: your logistic approach with interactions
  # If override: set all=1 or all=0
  logit_treat <- -1.2 +
    0.015*age +
    0.2*sex +
    0.4*cirrhosis +
    0.5*ckd +
    0.3*diabetes +
    0.2*hypertension +
    0.3*ment_ill +
    0.2*substance_abuse +
    0.4*(ckd*cirrhosis) +
    0.01*baseline_gfr

  # region effect
  logit_treat[region_raw=="S"] <- logit_treat[region_raw=="S"] - 0.2
  logit_treat[region_raw=="W"] <- logit_treat[region_raw=="W"] + 0.1

  # race effect
  for (i in seq_len(n)) {
    if (!is.na(race_raw[i])) {
      if (race_raw[i]=="black")    logit_treat[i] <- logit_treat[i] + 0.1
      if (race_raw[i]=="asian")    logit_treat[i] <- logit_treat[i] + 0.2
      if (race_raw[i]=="other")    logit_treat[i] <- logit_treat[i] - 0.1
    } else {
      logit_treat[i] <- logit_treat[i] + runif(1, -0.05, 0.05)
    }
  }

  if (treat_override == "random") {
    p_treat   <- plogis(logit_treat)
    treatment <- rbinom(n,1,p_treat)
  } else if (treat_override == "all_treated") {
    treatment <- rep(1, n)
  } else {
    treatment <- rep(0, n)
  }

  #### 3) Time-to-event with HR(treatment) ~ 1.06 ####
  # e.g. baseline hazard=0.02
  # We'll do main effect + small interactions, etc.
  # If you want same outcome model for all scenarios, do exactly that:
  # The only difference is we forcibly set 'treatment' above.
  lp_outcome <- -2.5 +
    0.03*age +
    0.6*ckd +
    0.4*cirrhosis +
    0.2*ment_ill +
    0.3*diabetes +
    0.1*substance_abuse +
    0.06*treatment +
    0.04*(treatment*ckd)

  rate_i <- 0.02*exp(lp_outcome)

  time_to_event <- rexp(n, rate=rate_i)
  censor_time   <- rexp(n, rate=1/100)
  observed_time <- pmin(time_to_event, censor_time)
  event         <- as.integer(time_to_event <= censor_time)

  #### 4) Build DF ####
  df <- data.frame(
    age=age, sex=sex, race=race_raw, region=region_raw,
    mental_illness=ment_ill, substance_abuse=substance_abuse, smoking=smoking,
    over_obesity=over_obesity, diabetes=diabetes, hypertension=hypertension,
    cirrhosis=cirrhosis, ckd=ckd, hiv=hiv, copd=copd, cancer=cancer,
    baseline_gfr=baseline_gfr, bmi=bmi,
    treatment=treatment,
    time=observed_time,
    event=event
  )

  #### 5) Some missingness (optional) ####
  set.seed(999)
  # region missing 5% except MW =>1%
  p_region_miss <- ifelse(df$region=="MW", 0.01, 0.05)
  reg_ix        <- runif(n) < p_region_miss
  df$region[reg_ix] <- NA

  # 10% missing ckd
  p_ckd_miss <- plogis(-3 + 1.0*df$ckd + 0.5*df$cirrhosis)
  ckd_ix     <- runif(n) < p_ckd_miss
  df$ckd[ckd_ix] <- NA

  return(df)
}


df0 <- simulate_data(n=5000, treat_override="all_untreated")
df1 <- simulate_data(n=5000, treat_override="all_treated")

# Combine them, label group=0 or 1
df0$group <- 0
df1$group <- 1
df_combined <- rbind(df0, df1)

# Now we have 10k rows, half are forced treatment=0, half=1
fit <- coxph(Surv(time, event) ~ group, data=df_combined)
summary(fit)  # This should recover the "true" HR ~ 1.06 (assuming big enough n)


# Using the function from above:

# 1) Everyone untreated
df0 <- simulate_data(n=5000, treat_override="all_untreated")

# 2) Everyone treated
df1 <- simulate_data(n=5000, treat_override="all_treated")

# Suppose we look at 1-year risk
t_star <- 365

risk_untreated <- mean(df0$time <= t_star & df0$event==1)
risk_treated   <- mean(df1$time <= t_star & df1$event==1)
risk_diff      <- risk_treated - risk_untreated
risk_ratio     <- risk_treated / risk_untreated

cat("1-year risk (untreated):", round(risk_untreated,4), "\n")
cat("1-year risk (treated):  ", round(risk_treated,4),   "\n")
cat("Risk difference:        ", round(risk_diff,4),      "\n")
cat("Risk ratio:        ", round(risk_diff,4),      "\n")

