# =============================================================
#  Simulated HCV‑AKI cohort – augmented to mirror the real Gilead data
#  Updated logit_treat to achieve ≈36 % SOF prevalence with good PS overlap
#  Source docs: rapid‑analysis spec & SAS dry‑run report citeturn6file0turn6file2
# =============================================================

# ---- 1. Libraries & set‑up --------------------------------------------
library(tidyverse)
library(missForest)
library(survival)
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

set.seed(12345)
h0 <- 0.02
lp_event <- with(cohort,
                 -2.8 + 0.03*age + 0.7*ckd + 0.5*cirrhosis + 0.3*heart_failure +
                   0.25*nsaid + 0.20*contrast)

#+ 0.7*treatment)
rate_i <- h0*exp(lp_event)

table(cohort$treatment)

haz_fun <- function(t, r, tr, HR_early = 1.5, HR_late = 0.7) {
  if (t <= 90) {                       # early segment
    r * ifelse(tr == 1, HR_early, 1)
  } else {                             # late segment
    r * ifelse(tr == 1, HR_late,  1)
  }
}

# ---------- simple piece-wise exponential RNG -----------------
rpexp1 <- function(n, rate1, rate2, tau = 90) {
  if (length(rate1) == 1) rate1 <- rep(rate1, n)
  if (length(rate2) == 1) rate2 <- rep(rate2, n)

  u  <- runif(n)
  p1 <- 1 - exp(-rate1 * tau)                 # P(T ≤ τ)

  out <- numeric(n)

  ## events in the first segment
  early        <- u <= p1
  out[early]   <- -log(1 - u[early]) / rate1[early]

  ## events after τ
  late         <- !early
  out[late]    <- tau -
    log((1 - u[late]) / (1 - p1[late])) / rate2[late]

  out
}

lambda1 <- haz_fun(30, rate_i, cohort$treatment)   # hazard 0–90 d
lambda2 <- haz_fun(91, rate_i, cohort$treatment)   # hazard >90 d

cohort$event_time <- rpexp1(nrow(cohort),
                            rate1 = lambda1,
                            rate2 = lambda2,
                            tau   = 90)


prop.table(table(cohort$event_time < 90))
summary(cohort$event_time )
summary(cohort$event_time[cohort$treatment==0])
summary(cohort$event_time[cohort$treatment==1])

cohort %>% group_by(treatment) %>%
  summarise(
    event_rate = mean(event_time),
    event_rate_90 = mean(event_time <= 90))

cohort$censor_admin <- pmin(rexp(nrow(cohort),1/100), max_follow)
cohort$switch       <- rbinom(nrow(cohort),1,.03)
cohort$censor_switch<- if_else(cohort$switch==1, cohort$tx_days + risk_window, max_follow)

cohort$follow_time <- pmin(cohort$event_time, cohort$censor_admin, cohort$censor_switch)
cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

#check hr
d90 <- cohort[cohort$follow_time<90,]
fit90 <- coxph(Surv(follow_time, event) ~ treatment, data = d90)
hr90  <- exp(coef(fit90))
hr90

fit_full <- coxph(Surv(follow_time, event) ~ treatment, data = cohort)
hr_full  <- exp(coef(fit_full))
hr_full


# ---- 10. Analysis dataset ---------------------------------------------
ana_vars <- cohort %>%
  dplyr::select(-enroll_days, -prior_aki, -prior_sof, -prior_nonsof,
                -tx_days, -event_time, -censor_admin, -censor_switch)

#Could evaluate missingness

# # ---- 11. Missingness (region 5 %, CKD 10 %) ---------------------------
# set.seed(123)
# ana_vars$region[ runif(nrow(ana_vars)) < 0.05 ] <- NA
# ana_vars$ckd[ runif(nrow(ana_vars)) < 0.10 ]    <- NA
#
# # ---- 12. Impute --------------------------------------------------------
# imp_vars <- c("age","race","region","ckd","cirrhosis","hiv","diabetes",
#               "hypertension","bmi")
# imp_in <- ana_vars[ , imp_vars ] %>% mutate(across(c(race,region), as.factor))
# imp_out <- missForest(as.data.frame(imp_in), verbose = TRUE)$ximp
# ana_vars[ , imp_vars] <- imp_out

# ---- 13. Save ----------------------------------------------------------
write_csv(ana_vars, "data/sim_hcv_aki_aug.csv")


fit_full <- coxph(Surv(follow_time, event) ~ treatment, data = ana_vars)
hr_full  <- exp(coef(fit_full))
hr_full



############################################################
##  Table-4 replica: basic rate parameters by exposure    ##
############################################################
library(tidyverse)
library(knitr)

## ---------- 1.  get a dataset (or skip if you already have `dat`) ------
## set.seed(123)
## dat <- generate_hcv_data(N = 70000,          # calibrated generator
##                          h0 = 3.5e-4, HR_early = 1.25,
##                          HR_late = 1.00, p_sof = 0.36)

## ---------- 2.  restrict follow-up to ≤180 d  --------------------------
tmax <- 180                               # analysis window (days)

dat180 <- ana_vars %>%
  mutate(fup_d  = pmin(follow_time, tmax),             # truncated time
         event180 = as.integer(event == 1 & follow_time <= tmax),
         fup_yrs = fup_d/365)

## ---------- 3.  arm-wise tallies --------------------------------------
tab <- dat180 %>%
  group_by(treatment) %>%
  summarise(
    `Number of patients`   = n(),
    `Person-years`         = sum(fup_yrs),
    `Number of events`     = sum(event180),
    .groups = "drop"
  ) %>%
  mutate(
    `Rate / 1,000 PY`      = `Number of events` / `Person-years` * 1e3
  )

## ---------- 4.  effect estimates  -------------------------------------
rate0 <- tab$`Rate / 1,000 PY`[tab$treatment==0]
rate1 <- tab$`Rate / 1,000 PY`[tab$treatment==1]
py0   <- tab$`Person-years`    [tab$treatment==0]
py1   <- tab$`Person-years`    [tab$treatment==1]
ev0   <- tab$`Number of events`[tab$treatment==0]
ev1   <- tab$`Number of events`[tab$treatment==1]

## Rate ratio  + Wald 95 % CI
RR      <- rate1 / rate0
se_logRR<- sqrt(1/ev1 + 1/ev0)
ciRR    <- exp(log(RR) + c(-1,1)*1.96*se_logRR)

## Rate difference  + Wald 95 % CI
RD      <- rate1 - rate0
se_RD   <- sqrt( ev1/(py1^2) + ev0/(py0^2) ) * 1e3   # converted to /1 000 PY
ciRD    <- RD + c(-1,1)*1.96*se_RD

## ---------- 5.  assemble final table -----------------------------------
out <- tab %>%
  mutate(
    treatment = factor(treatment, levels = c(0,1),
                       labels = c("non-SOF DAAs", "SOF-DAAs"))
  ) %>%
  dplyr::select(Parameter = treatment,
                `Number of patients`,
                `Person-years`,
                `Number of events`,
                `Rate / 1,000 PY`) %>%
  mutate(`Number of patients`=as.character(`Number of patients`),
         `Person-years`      = sprintf("%.1f", `Person-years`),
         `Number of events`  = as.character(`Number of events`),
         `Rate / 1,000 PY`   = sprintf("%.2f", `Rate / 1,000 PY`)) %>%
  bind_rows(
    tibble(Parameter = "Rate ratio (vs. referent; 95% CI)",
           `Number of patients` = "Referent",
           `Person-years`       = "",
           `Number of events`   = "",
           `Rate / 1,000 PY`    = sprintf("%.2f (%.2f, %.2f)", RR, ciRR[1], ciRR[2])),
    tibble(Parameter = "Rate difference / 1,000 PY (95% CI)",
           `Number of patients` = "Referent",
           `Person-years`       = "",
           `Number of events`   = "",
           `Rate / 1,000 PY`    = sprintf("%.2f (%.2f, %.2f)", RD, ciRD[1], ciRD[2])))


out
