# ─────────────────────────────────────────────────────────────────────────────
# Parallel sim: LMTP SDR vs PS-matched Cox HR
# ─────────────────────────────────────────────────────────────────────────────
library(data.table)
library(dplyr)
library(tidyr)
library(lmtp)
library(SuperLearner)
library(furrr)
library(here)
library(MatchIt)
library(survival)
library(broom)

# 1) Parallel setup
future::plan(multisession, workers = 4)     # Windows-friendly

# 2) Simulation parameters
n_iter      <- 50       # number of replications
sample_size <- 1000     # sample per iteration
tau         <- 180      # follow-up window

# 3) Load full data once
fulldat <- fread(here("data/sim_hcv_aki_complex.csv"))

# 4) PS formula for matching
ps_formula <- treatment ~ age + sex_male + race + region +
  cirrhosis + hiv + ckd +
  diabetes + hypertension + sepsis + bmi

# 5) Single-iteration function
run_one <- function(i) {
  set.seed(1234 + i)

  # 5a) Sample
  dat <- fulldat %>% sample_n(sample_size, replace = TRUE)

  # 5b) Propensity-score matching & Cox HR
  ps_mod <- glm(ps_formula, data=dat, family=binomial)
  dat$ps <- predict(ps_mod, type="response")

  m_out <- matchit(ps_formula, data=dat, method="nearest",
                   caliper=0.2, ratio=1)
  m_dat <- match.data(m_out)

  cox_fit <- coxph(Surv(follow_time, event) ~ treatment, data=m_dat)
  hr_tab  <- tidy(cox_fit, exponentiate=TRUE, conf.int=TRUE) %>%
    filter(term=="treatment") %>%
    select(estimate, conf.low, conf.high, p.value)

  # 5c) Save
  saveRDS(
    list(
      iter       = i,
      cox_hr     = hr_tab$estimate,
      cox_ci_low = hr_tab$conf.low,
      cox_ci_high = hr_tab$conf.high,
      pval = hr_tab$p.value
    ),
    file = here("results","sim_results", paste0("ps_", i, ".rds"))
  )

  return(data.frame(
    iter       = i,
    cox_hr     = hr_tab$estimate,
    cox_ci_low = hr_tab$conf.low,
    cox_ci_high= hr_tab$conf.high,
    pval = hr_tab$p.value
  ))
}

#DEBUG!

# # 6) Run sims in parallel
# res_df <- future_map_dfr(1:n_iter, run_one, .progress=TRUE)

res_df = NULL
for(i in 1:n_iter){
  res = run_one(i)
  res_df = bind_rows(res_df, res)
}


# 7) Inspect aggregate results
print(head(res_df))
write.csv(res_df, here("results","summary_ps.csv"), row.names=FALSE)
