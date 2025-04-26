# =============================================================
#  Simulated HCV‑AKI cohort – augmented to mirror the real Gilead data
#  Updated logit_treat to achieve ≈36 % SOF prevalence with good PS overlap
#  Source docs: rapid‑analysis spec & SAS dry‑run report citeturn6file0turn6file2
# =============================================================

# ---- 1. Libraries & set‑up --------------------------------------------
library(tidyverse)
library(missForest)
set.seed(42)

# ---- 2. Parameters -----------------------------------------------------
N            <- 20000      # base N (allows exclusions)
max_follow   <- 180        # 6‑month risk window
risk_window  <- 30         # 30‑d grace after EOT / switch

# ---- 3. Demography & enrolment ----------------------------------------
raw <- tibble(
  id       = 1:N,
  age      = pmax(rnorm(N, 48, 13), 18),
  sex_male = rbinom(N, 1, 0.58),
  race     = sample(c("white","black","hispanic","asian","other"), N, TRUE,
                    prob = c(.48,.14,.06,.02,.30)),
  region   = sample(c("NE","MW","S","W"), N, TRUE, prob = c(.20,.18,.37,.25)),
  enroll_days = rpois(N, 420)                    # prior continuous coverage
) %>%
  mutate(race = replace(race, runif(N) < .20, NA))      # 20 % missing race

# ---- 4. Clinical history & behaviour ----------------------------------
raw <- raw %>%
  mutate(
    ckd     = rbinom(N,1,.08),
    prior_aki = rbinom(N,1,.05),
    heart_failure = rbinom(N,1,.07),
    sepsis  = rbinom(N,1,.03),
    dehydration = rbinom(N,1,.06),
    obstruction = rbinom(N,1,.04),
    cirrhosis = rbinom(N,1,.18),
    portal_htn = rbinom(N,1,.04),
    esld   = rbinom(N,1,.02),
    hiv    = rbinom(N,1,.04),
    diabetes = rbinom(N,1,.20),
    hypertension = rbinom(N,1,.45),
    bmi    = rnorm(N,28,5),
    overweight_obese = rbinom(N,1,.20),
    smoking = rbinom(N,1,.40),
    alcohol = rbinom(N,1,.18),
    substance_abuse = rbinom(N,1,.25),
    cancer = rbinom(N,1,.08),
    chemo  = rbinom(N,1,.01)
  )

# ---- 5. Medications ----------------------------------------------------
raw <- raw %>%
  mutate(
    nsaid   = rbinom(N,1,.25),
    acearb  = rbinom(N,1,.30),
    diuretic = rbinom(N,1,.22),
    aminoglycoside = rbinom(N,1,.05),
    contrast = rbinom(N,1,.08),
    statin  = rbinom(N,1,.15),
    aspirin = rbinom(N,1,.10),
    beta_blocker = rbinom(N,1,.14),
    ccb     = rbinom(N,1,.16),
    art     = rbinom(N,1,.05)
  )

# ---- 6. Prior DAA exposure (wash‑out) ----------------------------------
raw <- raw %>%
  mutate(
    prior_sof    = rbinom(N,1,.05),
    prior_nonsof = rbinom(N,1,.05)
  )

# ---- 7. Apply baseline exclusions -------------------------------------
cohort <- raw %>%
  filter(
    enroll_days >= 365,
    age >= 18,
    prior_aki == 0,
    !(prior_sof == 1 | prior_nonsof == 1)
  )

# ---- 8. Treatment assignment (new!, tuned to real data) ---------------
# Goal: marginal P(SOF)=0.36, moderate channeling, good PS overlap

lp <- with(cohort,
           0.015*age +
             0.30*cirrhosis + 0.35*portal_htn + 0.25*ckd +
             0.20*hiv + 0.15*substance_abuse + 0.10*diabetes + 0.05*hypertension +
             -0.10*cancer + -0.05*overweight_obese +
             if_else(region=="W", 0.05, if_else(region=="S", -0.05, 0)) +
             case_when(race=="black" ~ 0.05,
                       race=="asian" ~ 0.10,
                       race=="other" ~ -0.05,
                       TRUE ~ 0) +
             rnorm(nrow(cohort), 0, 0.6)            # unmeasured drivers
)

alpha0 <- qlogis(0.36) - mean(lp)        # calibrate to 36 % prevalence
logit_ps <- alpha0 + lp
p_treat  <- plogis(logit_ps)
# enforce overlap to avoid extreme weights
p_treat  <- pmin(pmax(p_treat, 0.05), 0.95)
cohort$treatment <- rbinom(nrow(cohort), 1, p_treat)   # 1 = SOF

# ---- 9. Treatment duration & follow‑up --------------------------------
cohort$tx_days <- if_else(cohort$treatment==1,
                          rpois(nrow(cohort), 84),
                          rpois(nrow(cohort), 70))

h0 <- 0.02
lp_event <- with(cohort,
                 -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis + 0.3*heart_failure +
                   0.25*nsaid + 0.20*contrast + 0.7*treatment)
rate_i <- h0*exp(lp_event)

haz_fun <- function(t,r,tr){ ifelse(t<=90, r*if_else(tr==1,1.5,1), r*if_else(tr==1,0.7,1)) }

u <- runif(nrow(cohort))
cohort$event_time <- if_else(
  u < 1 - exp(-haz_fun(90,rate_i,cohort$treatment)*90),
  -log(1-u)/haz_fun(90,rate_i,cohort$treatment),
  90 + rexp(nrow(cohort), rate = haz_fun(91,rate_i,cohort$treatment))
)

cohort$censor_admin <- pmin(rexp(nrow(cohort),1/100), max_follow)
cohort$switch       <- rbinom(nrow(cohort),1,.03)
cohort$censor_switch<- if_else(cohort$switch==1, cohort$tx_days + risk_window, max_follow)

cohort$follow_time <- pmin(cohort$event_time, cohort$censor_admin, cohort$censor_switch)
cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

# ---- 10. Analysis dataset ---------------------------------------------
ana_vars <- cohort %>%
  select(-enroll_days, -prior_aki, -prior_sof, -prior_nonsof,
         -tx_days, -event_time, -censor_admin, -censor_switch)

#Could evaluate missingness

# ---- 11. Missingness (region 5 %, CKD 10 %) ---------------------------
set.seed(123)
ana_vars$region[ runif(nrow(ana_vars)) < 0.05 ] <- NA
ana_vars$ckd[ runif(nrow(ana_vars)) < 0.10 ]    <- NA

# ---- 12. Impute --------------------------------------------------------
imp_vars <- c("age","race","region","ckd","cirrhosis","hiv","diabetes",
              "hypertension","bmi")
imp_in <- ana_vars[ , imp_vars ] %>% mutate(across(c(race,region), as.factor))
imp_out <- missForest(as.data.frame(imp_in), verbose = TRUE)$ximp
ana_vars[ , imp_vars] <- imp_out

# ---- 13. Save ----------------------------------------------------------
write_csv(ana_vars, "data/sim_hcv_aki_aug.csv")

# Data ready:  `ana_vars`
