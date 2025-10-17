# Load required packages
library(data.table)
library(dplyr)
library(tidyr)
library(lmtp)
library(SuperLearner)
options(future.globals.maxSize = 6 * 1024^3)

# 1. Load data
fulldat <- read.csv(here::here("data/sim_hcv_aki_complex.csv"))

start_time <- proc.time()


#subset for speed temporarily
dat <- fulldat %>% sample_n(1000, replace=TRUE)

# 2. Define follow-up horizon
tau <- 180  # e.g., 180-day cumulative incidence window

# 3. Derive time-to-event structure
dat <- dat %>%
  mutate(
    aki_event = as.integer(event == 1 & follow_time <= tau),
    time_to_aki = if_else(aki_event == 1, follow_time, NA_real_),
    cens_event = as.integer(event == 0 & follow_time < tau),
    time_to_censor = if_else(cens_event == 1, follow_time, NA_real_)
  )

#sanity checks
summary(dat$time_to_aki)
table(dat$aki_event, dat$cens_event)
prop.table(table(dat$treatment))


# 4. Expand to wide-format Y and C
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

# 5. Rename Y and C columns
Y_cols <- paste0("Yvec", 1:tau)
C_cols <- paste0("Cvec", 1:tau)
setnames(wide, old = paste0("Y_vec", 1:tau), new = Y_cols)
setnames(wide, old = paste0("C_vec", 1:tau), new = C_cols)

# 6. Carry forward AKI once it happens
wide <- event_locf(wide, outcomes = Y_cols)

wide = as.data.frame(wide)



# 7. Run LMTP SDR
res <- lmtp_sdr(
  data = wide,
  trt = "treatment",           # exposure variable
  outcome = Y_cols,
  cens = C_cols,
  baseline = c("age", "sex_male", "ckd", "diabetes", "hypertension"),  # minimal covariate set; customize as needed
  shift = static_binary_on,   # set A=1 (treat all with SOF)
  outcome_type = "survival",
  folds = 2, # for speed
  learners_trt = c("SL.glm"),
  learners_outcome = c("SL.glm")
)

print(res)

res_0 <- lmtp_sdr(
  data = wide,
  trt = "treatment",           # exposure variable
  outcome = Y_cols,
  cens = C_cols,
  baseline = c("age", "sex_male", "ckd", "diabetes", "hypertension"),  # minimal covariate set; customize as needed
  shift = static_binary_off,   # set A=1 (treat all with SOF)
  outcome_type = "survival",
  folds = 2, # for speed
  learners_trt = c("SL.glm"),
  learners_outcome = c("SL.glm")
)

print(res_0)

# # 8. Estimate contrast: treat-all vs. observed
# res_nat <- lmtp_sdr(
#   data = wide,
#   trt = "treatment",
#   outcome = Y_cols,
#   cens = C_cols,
#   baseline = c("age", "sex_male", "ckd", "diabetes", "hypertension"),
#   shift = NULL,               # observed treatment
#   outcome_type = "survival",
#   folds = 5,
#   learners_trt = c("SL.glm"),
#   learners_outcome = c("SL.glm")
# )

# 9. Estimate cumulative incidence ratio
# lmtp_contrast(res, ref = res_nat, type = "rr")


lmtp_contrast(res, ref = res_0, type = "rr")


end_time <- proc.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

