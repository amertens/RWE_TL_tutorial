
library(survival)
library(tmle)
library(survival)
library(survey)
library(here)


source(here("DGP.R"))   # assumes the updated function is in your working dir

make_landmark <- function(dat, t0){
  dat |> dplyr::mutate(
    Y     = dplyr::case_when(event == 1 & follow_time <= t0 ~ 1,
                             follow_time  >= t0             ~ 0,
                             TRUE ~ NA_real_),
    Delta = ifelse(is.na(Y), 0, 1)
  )
}




## 4.2 Single-dataset analysis


N = 2000
horizons = c(90)
analyze_once <- function(N = 20000, horizons = c(30, 90, 180)){
  dat   <- generate_hcv_data(N = N, np_hazard    = TRUE, dep_censor    = TRUE, complexity    = TRUE)
  Wvars <- c("age","sex_male","ckd","cirrhosis","diabetes","hypertension",
             "bmi","race","region")

  lapply(horizons, function(t0){
    d0 <- make_landmark(dat, t0)
    W  <- d0[, Wvars]

    fit <- tmle(Y = d0$Y, A = d0$treatment, W = W, Delta = d0$Delta,
                family = "binomial",
                Q.SL.library = c("SL.glm","SL.mean"),
                g.SL.library = "SL.glm",
                g.Delta.SL.library = "SL.glm")
    R1 <- fit$estimates$EY1$psi
    R0 <- fit$estimates$EY0$psi

    RD <- R1 - R0
    RR <- R1 / R0

    hr_crude <- exp(coef(coxph(Surv(follow_time,event) ~ treatment, data = dat)))

    ps  <- glm(treatment ~ ., data = dat[, c("treatment", Wvars)], family = binomial)$fitted
    w   <- ifelse(dat$treatment==1, 1/ps, 1/(1-ps))
    des <- svydesign(ids = ~1, weights = ~w, data = dat)
    hr_iptw <- exp(coef(svycoxph(Surv(follow_time,event) ~ treatment, design = des)))

    ps_form <- as.formula(paste("treatment ~", paste(Wvars, collapse = "+")))

    m.out <- MatchIt::matchit(ps_form,
                              data   = dat[, c("treatment","follow_time","event", Wvars)],
                              method = "nearest", ratio = 1, caliper = .20)

    d.m   <- MatchIt::match.data(m.out, data = "all")

    hr_psm <- exp(coef(coxph(Surv(follow_time, event) ~ treatment,
                             data = d.m, cluster = subclass)))

    tibble(t0 = t0, R1, R0, RD, RR, hr_crude, hr_iptw, hr_psm)
  }) |> dplyr::bind_rows()
}

set.seed(2025)
#res=analyze_once(N = 20000, horizons = c(90))

B <- 10   # iterations; adjust as needed
res <- purrr::map_dfr(1:B, ~analyze_once(N = 4000, horizons = c(90)))

sim_sum <- res |>
  group_by(t0) |>
  summarise(across(R1:hr_iptw, list(mean = mean, sd = sd), .names = "{col}_{fn}"))
knitr::kable(sim_sum, digits = 3)



# ---------- performance summary ----------------------------------------
perf_stats <- function(res, truth,
                       est_cols   = c("RD", "RR", "hr_crude", "hr_psm"),
                       lwr_suffix = "_lwr", upr_suffix = "_upr"){
  truth_vec <- truth[est_cols]


  ## helper that returns a one-row tibble for a single estimand
  one_stat <- function(col) {
    x <- res[[col]]
    tibble(
      estimand = col,
      mean     = mean(x),
      sd       = sd(x),                    # oracle SE
      bias     = mean(x) - truth[[col]],
      cover    = {
        lwr <- res[[paste0(col, lwr_suffix)]]
        upr <- res[[paste0(col, upr_suffix)]]
        if (!is.null(lwr) && !is.null(upr))
          mean(lwr <= truth[[col]] & upr >= truth[[col]]) * 100
        else NA_real_
      }
    )
  }

  purrr::map_dfr(est_cols, one_stat)

}





truth <- readRDS(file = here("results/truth.rds"))

### --- 2. pick estimand columns ----------------------------------------
estimands <- intersect(names(res),
                       c("RD","RR","hr_crude","hr_psm",
                         "hr_iptw","hr_tmle"))

truth <- truth %>% mutate(hr_crude = HR,
                          hr_psm = HR,
                          hr_iptw = HR)



# ---------- performance summary ----------------------------------------
perf_stats <- function(res, truth,
                       est_cols   = c("RD", "RR", "hr_crude", "hr_psm"),
                       lwr_suffix = "_lwr", upr_suffix = "_upr"){
  truth_vec <- truth[est_cols]


  ## helper that returns a one-row tibble for a single estimand
  one_stat <- function(col) {
    x <- res[[col]]
    tibble(
      estimand = col,
      mean     = mean(x),
      sd       = sd(x),                    # oracle SE
      bias     = mean(x) - truth[[col]],
      cover    = {
        lwr <- res[[paste0(col, lwr_suffix)]]
        upr <- res[[paste0(col, upr_suffix)]]
        if (!is.null(lwr) && !is.null(upr))
          mean(lwr <= truth[[col]] & upr >= truth[[col]]) * 100
        else NA_real_
      }
    )
  }

  purrr::map_dfr(est_cols, one_stat)

}

# --------------------  use ---------------------------------------------
stats_out <- perf_stats(res  = res,          # your replicate tibble
                        truth = truth)
print(stats_out, digits = 3)
