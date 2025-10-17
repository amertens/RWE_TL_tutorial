# Load required packages
library(data.table)
library(dplyr)
library(tidyr)
library(lmtp)
library(SuperLearner)
library(furrr)
library(here)

# Setup parallel plan for Windows
future::plan(multisession, workers = 20)  # adjust workers to your CPU

# Define simulation parameters
n_iter <- 200            # number of iterations
sample_size <- 10000     # sample size per iteration
tau <- 180              # max follow-up

# Load full data once
fulldat <- fread(here("data/sim_hcv_aki_complex.csv"))

# Main wrapper function for one iteration
run_iteration <- function(i) {
  set.seed(1234 + i)
  dat <- fulldat %>% sample_n(sample_size, replace = TRUE)

  dat <- dat %>%
    mutate(
      aki_event = as.integer(event == 1 & follow_time <= tau),
      time_to_aki = if_else(aki_event == 1, follow_time, NA_real_),
      cens_event = as.integer(event == 0 & follow_time < tau),
      time_to_censor = if_else(cens_event == 1, follow_time, NA_real_)
    )

  make_Y <- function(t, e) as.integer(1:tau >= t & e == 1)
  make_C <- function(t, e) {
    v <- rep(1, tau)
    if (e == 1) v[(t + 1):tau] <- 0
    v
  }

  wide <- dat %>%
    rowwise() %>%
    mutate(
      Y_vec = list(make_Y(time_to_aki, aki_event)),
      C_vec = list(make_C(time_to_censor, cens_event))
    ) %>%
    unnest_wider(Y_vec, names_sep = "") %>%
    unnest_wider(C_vec, names_sep = "")

  Y_cols <- paste0("Yvec", 1:tau)
  C_cols <- paste0("Cvec", 1:tau)
  setnames(wide, paste0("Y_vec", 1:tau), Y_cols)
  setnames(wide, paste0("C_vec", 1:tau), C_cols)

  wide <- event_locf(wide, outcomes = Y_cols)
  wide <- as.data.frame(wide)

  baseline_vars <- c("age", "sex_male", "ckd", "diabetes", "hypertension")

  res1 <- lmtp_sdr(
    data = wide,
    trt = "treatment",
    outcome = Y_cols,
    cens = C_cols,
    baseline = baseline_vars,
    shift = static_binary_on,
    outcome_type = "survival",
    folds = 2,
    learners_trt = c("SL.glm"),
    learners_outcome = c("SL.glm")
  )

  res0 <- lmtp_sdr(
    data = wide,
    trt = "treatment",
    outcome = Y_cols,
    cens = C_cols,
    baseline = baseline_vars,
    shift = static_binary_off,
    outcome_type = "survival",
    folds = 2,
    learners_trt = c("SL.glm"),
    learners_outcome = c("SL.glm")
  )

  contrast <- lmtp_contrast(res1, ref = res0, type = "rr")

  saveRDS(list(res_on = res1, res_off = res0, contrast = contrast),
          file = here("results","sim_results", paste0("lmtp_sim_", i, ".rds")))

  return(contrast)
}


# Run parallel loop
start_time <- proc.time()
all_contrasts <- future_map(1:n_iter, run_iteration, .progress = TRUE)
end_time <- proc.time()
print(end_time - start_time)

saveRDS(all_contrasts, file = here("results","sim_results", "lmtp_sim_all.rds"))



