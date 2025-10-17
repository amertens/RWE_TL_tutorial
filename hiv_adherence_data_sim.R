suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(lubridate); library(data.table)
})
set.seed(1)

n_pat = 3000
start_date = as.Date("2019-01-01")
end_date   = as.Date("2025-12-31")
recruit_start = as.Date("2020-01-01")
recruit_end   = as.Date("2025-01-31")
max_follow_days = 365*4   # allow up to 4y per person in raw data
p_lab_link = 0.85


sim_claims_like <- function(n_pat = 3000,
                            start_date = as.Date("2019-01-01"),
                            end_date   = as.Date("2025-12-31"),
                            recruit_start = as.Date("2020-01-01"),
                            recruit_end   = as.Date("2025-01-31"),
                            max_follow_days = 365*4,   # allow up to 4y per person in raw data
                            p_lab_link = 0.85) {

  # --- patients ---
  pt <- tibble(
    pat_id = sprintf("P%06d", 1:n_pat),
    sex = rbinom(n_pat, 1, 0.55),                 # 1=male
    dob = as.Date("1955-01-01") + sample(0:20000, n_pat, TRUE),
    lab_link = rbinom(n_pat, 1, p_lab_link)
  ) %>%
    mutate(age_2020 = as.integer(time_length(interval(dob, as.Date("2020-01-01")), "years")))

  # --- eligibility (enrollment) segments ---
  # each person: 1–3 segments; some gaps; lengths 6–36 months
  segs_per <- sample(1:3, n_pat, TRUE, prob = c(0.7, 0.25, 0.05))
  elig <- map2_dfr(pt$pat_id, segs_per, function(id, k){
    s <- sort(sample(seq(start_date, end_date - 180, by="1 month"), k))
    len <- sample(180:1080, k, TRUE)  # 6–36 months
    tibble(pat_id = id,
           enr_start = s,
           enr_end   = pmin(s + len, end_date))
  }) %>%
    arrange(pat_id, enr_start)

  # --- baseline characteristics for disease process ---
  base <- pt %>%
    transmute(pat_id, sex, age0 = pmin(85, pmax(18, rnorm(n(), 42, 10))),
              cd4_0 = pmax(20, rnorm(n(), 500, 180)),
              vl0   = exp(rnorm(n(), log(1e5), 0.7))) %>%
    mutate(eta_adher = 0.5*scale(cd4_0)[,1] - 0.2*(log10(vl0)-5) + 0.2*sex + rnorm(n(),0,0.7))

  # --- choose an index window with >=12m prior eligibility ---
  # find first time a person has an eligibility run >= 24 months to host baseline+follow-up
  elig_long <- elig %>%
    group_by(pat_id) %>%
    mutate(run_len = as.integer(enr_end - enr_start)) %>%
    ungroup() %>%
    filter(run_len >= 365*2) %>%
    mutate(idx_window_start = pmax(enr_start + 365, recruit_start),
           idx_window_end   = pmin(enr_end - 365, recruit_end)) %>%   # allow 12m baseline and 12m follow-up at least
    filter(idx_window_start <= idx_window_end)

  # only keep patients with a valid index window
  keep_ids <- unique(elig_long$pat_id)
  base <- base %>% filter(pat_id %in% keep_ids)
  pt   <- pt %>% filter(pat_id %in% keep_ids)

  # --- simulate index regimen & date ---
  # regimen A vs B probability depends on sex, age, baseline adherence propensity
  idx_df <- elig_long %>%
    group_by(pat_id) %>%
    slice_sample(n = 1) %>% ungroup() %>%
    left_join(base, by="pat_id") %>%
    mutate(pA = plogis(-0.1 + 0.25*sex - 0.1*(age0-45)/10 + 0.4*eta_adher),
           regimen = rbinom(n(), 1, pA),   # 1=A, 0=B
           index_date = as.Date(runif(n(), as.numeric(idx_window_start), as.numeric(idx_window_end)), origin = "1970-01-01"))

  # --- pharmacy fills ---
  # new-user: **no ART fills in the prior 365 days** by construction
  gen_fills_one <- function(id, idx, reg, eta, enr_start, enr_end) {
    # adhere intensity -> gap distribution
    max_t <- min(idx + max_follow_days, enr_end)
    if (idx > max_t) return(NULL)
    # simulate primary regimen fills
    days_supply <- 30L
    fill_dates <- as.Date(idx)
    # daily hazard of refill = function of adherence propensity
    while (TRUE) {
      last <- tail(fill_dates, 1)
      # "gap" days: geometric-ish with mean depending on eta (lower gap if high adherence)
      mean_gap <- 10 - 4*plogis(eta)   # mean ~6–10
      gap <- pmax(0L, rpois(1, lambda = mean_gap))
      next_date <- last + days_supply + gap
      if (next_date > max_t) break
      fill_dates <- c(fill_dates, next_date)
      if (length(fill_dates) > 30) break
    }
    # possible switch (~8–12% by 12m), earlier if low adherence and high VL
    switch_flag <- rbinom(1, 1, prob = 0.12 - 0.08*plogis(eta)) == 1
    switch_date <- if (switch_flag) as.Date(idx + sample(90:360,1)) else as.Date(NA)
    # build rows
    tibble(
      pat_id = id,
      ndc = ifelse(reg==1, "NDC_A", "NDC_B"),
      drug_class = ifelse(reg==1, "ART_A", "ART_B"),
      fill_date = fill_dates,
      days_supply = 30L,
      quantity = 30L,
      index_regimen = reg,
      switch_date = switch_date
    )
  }

  rx <- pmap_dfr(idx_df %>% select(pat_id, index_date, regimen) %>%
                   left_join(elig %>% group_by(pat_id) %>% slice_min(enr_start, with_ties = FALSE) %>% ungroup(),
                             by="pat_id") %>%
                   left_join(elig %>% group_by(pat_id) %>% slice_max(enr_end, with_ties = FALSE) %>% ungroup(),
                             by="pat_id") %>%
                   left_join(base %>% select(pat_id, eta_adher), by="pat_id"),
                 ~ gen_fills_one(..1, ..2, ..3, ..6, ..4, ..5)) %>%
    arrange(pat_id, fill_date)

  # add post-switch fills to opposite class (simple pattern)
  rx_sw <- rx %>%
    group_by(pat_id) %>%
    summarise(switch_date = first(na.omit(switch_date)), index_regimen = first(index_regimen), .groups="drop") %>%
    filter(!is.na(switch_date)) %>%
    mutate(ndc = ifelse(index_regimen==1, "NDC_B", "NDC_A"),
           drug_class = ifelse(index_regimen==1, "ART_B", "ART_A")) %>%
    group_by(pat_id) %>%
    do({
      id <- .$pat_id[1]; sdt <- .$switch_date[1]; ndc <- .$ndc[1]; dc <- .$drug_class[1]
      nfill <- sample(2:6,1); days_supply <- 30L
      tibble(pat_id=id, ndc=ndc, drug_class=dc,
             fill_date = sdt + cumsum(sample(25:40, nfill, TRUE)),
             days_supply = 30L, quantity = 30L)
    }) %>% ungroup()

  rx <- bind_rows(rx, rx_sw) %>% arrange(pat_id, fill_date)

  # --- diagnoses (HIV codes) ---
  dx <- idx_df %>%
    mutate(n_hiv_dx = sample(1:3, n(), TRUE, c(0.1, 0.7, 0.2))) %>%
    rowwise() %>%
    do({
      tibble(pat_id = .$pat_id,
             svc_date = sort(.$index_date - sample(30:400, .$n_hiv_dx, TRUE)),
             icd10 = sample(c("B20","Z21"), .$n_hiv_dx, TRUE))
    }) %>% ungroup()

  # --- labs (VL + CD4); informative monitoring ---
  gen_labs_one <- function(id, idx, lab_on, base_cd4, base_vl, eta, regimen) {
    # baseline windows
    base_dates <- sort(idx - sample(0:180, size = sample(1:2,1), replace=TRUE))
    # follow-up: roughly q6 months, but probability depends on adherence & prior VL
    follow_dates <- idx + cumsum(sample(120:240, size = sample(2:6,1), replace=TRUE))
    all_dates <- unique(sort(c(base_dates, follow_dates)))
    if (!lab_on) all_dates <- base_dates[1]  # lab_link==0 → almost no labs

    # generate values
    log10_vl <- log10(base_vl); cd4 <- base_cd4
    out <- map_dfr(all_dates, function(d){
      # regimen at date: assume index regimen until switch (approx)
      # simple dynamic: adherence propensity lowers VL over time; poor adherence raises VL
      log10_vl <<- pmax(2.0, log10_vl + 0.05 - 0.25*regimen - 0.15*plogis(eta) + rnorm(1,0,0.25))
      cd4 <<- pmax(20, cd4 + rnorm(1, 5*plogis(eta), 40))
      tibble(pat_id=id, lab_date=d, loinc=c("HIVRNA","CD4"),
             value = c(10^rnorm(1, log10_vl, 0.2), cd4))
    })
    out
  }

  labs <- pmap_dfr(idx_df %>% select(pat_id, index_date, regimen) %>%
                     left_join(pt %>% select(pat_id, lab_link), by="pat_id") %>%
                     left_join(base %>% select(pat_id, cd4_0, vl0, eta_adher), by="pat_id"),
                   ~ gen_labs_one(..1, ..2, ..4==1, ..5, ..6, ..7, ..3)) %>%
    mutate(test = ifelse(loinc=="HIVRNA", "VL", "CD4"),
           result = ifelse(test=="VL", value, value)) %>%
    select(pat_id, lab_date, test, result)

  # --- death dates (low rate; older age higher risk) ---
  death <- idx_df %>%
    transmute(pat_id,
              death = rbinom(n(), 1, plogis(-6 + 0.03*(age0-50)))==1,
              death_date = ifelse(death, as.Date(index_date + sample(120:900,1)), NA)) %>%
    filter(death) %>% select(pat_id, death_date)

  list(pt=pt, elig=elig, rx=rx, labs=labs, dx=dx, death=death, idx_df=idx_df)
}

raw <- sim_claims_like()
str(raw, max.level = 1)



# helper: continuous enrollment check for a window
has_cont_enroll <- function(elig_tbl, id, start_d, end_d) {
  e <- elig_tbl %>% filter(pat_id==id) %>% arrange(enr_start)
  # coverage if any segment fully spans [start_d, end_d]
  any(e$enr_start <= start_d & e$enr_end >= end_d)
}

# classify ART class A/B fills (already coded in simulation)
rx_ab <- raw$rx %>% filter(drug_class %in% c("ART_A","ART_B"))

# first qualifying ART fill as **index**
first_art <- rx_ab %>% group_by(pat_id) %>% slice_min(fill_date, with_ties=FALSE) %>%
  ungroup() %>% rename(index_date = fill_date, index_class = drug_class) %>%
  mutate(A0 = as.integer(index_class=="ART_A"))

# baseline *no prior ART* in 365d
no_prior_art <- rx_ab %>%
  inner_join(first_art %>% select(pat_id, index_date), by="pat_id") %>%
  filter(fill_date < index_date & fill_date >= index_date - 365) %>%
  distinct(pat_id) %>% mutate(prior_art=1L)

# baseline labs in [-180,+30] around index (VL and CD4)
base_labs <- raw$labs %>%
  semi_join(first_art, by="pat_id") %>%
  inner_join(first_art %>% select(pat_id, index_date), by="pat_id") %>%
  filter(lab_date >= index_date - 180 & lab_date <= index_date + 30) %>%
  group_by(pat_id) %>%
  summarise(has_vl = any(test=="VL"),
            has_cd4= any(test=="CD4"), .groups="drop")

# lab link & at least one post-index VL
post_vl <- raw$labs %>%
  inner_join(first_art %>% select(pat_id, index_date), by="pat_id") %>%
  filter(lab_date > index_date) %>%
  group_by(pat_id) %>% summarise(has_post_vl = any(test=="VL"), .groups="drop")

# death date
death_dt <- raw$death

# apply inclusion/exclusion
cohort <- first_art %>%
  left_join(raw$pt %>% select(pat_id, sex, dob, lab_link), by="pat_id") %>%
  mutate(age = as.integer(time_length(interval(dob, index_date), "years"))) %>%
  left_join(base_labs, by="pat_id") %>%
  left_join(post_vl, by="pat_id") %>%
  left_join(death_dt, by="pat_id") %>%
  anti_join(no_prior_art, by="pat_id") %>%
  filter(age >= 18,
         lab_link==1,
         has_vl, has_cd4,
         has_post_vl) %>%
  # continuous enrollment: 365d baseline + desired follow-up horizon (here 12m)
  rowwise() %>%
  mutate(baseline_ok = has_cont_enroll(raw$elig, pat_id, index_date-365, index_date-1),
         follow_ok   = has_cont_enroll(raw$elig, pat_id, index_date,   index_date + 365)) %>%
  ungroup() %>% filter(baseline_ok, follow_ok)

nrow(cohort)


# block grid per person
make_blocks <- function(index_date, K=4, len=90) {
  tibble(k = 1:K,
         start = index_date + (k-1)*len,
         end   = index_date + k*len - 1L)
}

# compute covered days set from fills with carryover policy
covered_days <- function(fills, carry=c("full","cap30"), cap=30L) {
  carry <- match.arg(carry)
  if (nrow(fills)==0) return(as.Date(character()))
  fills <- fills %>% arrange(fill_date)
  # build intervals with or without capping
  if (carry=="full") {
    iv <- IRanges::IRanges(start = as.integer(fills$fill_date),
                           width = fills$days_supply)
    iv <- IRanges::reduce(iv)
    d <- as.Date(IRanges::as.data.frame(iv) %>%
                   rowwise() %>%
                   do(tibble(day = as.integer(.$start):as.integer(.$end))) %>%
                   pull(day), origin="1970-01-01")
    return(d)
  } else {
    # daily inventory with 30-day cap on stock
    min_day <- min(fills$fill_date)
    max_day <- max(fills$fill_date + fills$days_supply + 120)
    days <- seq(min_day, max_day, by="day")
    stock <- integer(length(days))
    date_idx <- match(fills$fill_date, days)
    for (i in seq_along(days)) {
      # add refills arriving today
      add <- sum(date_idx == i) * 1L
      if (add>0) stock[i] <- min(cap, (if (i>1) stock[i-1] else 0) + sum(fills$days_supply[date_idx==i]))
      else stock[i] <- if (i>1) stock[i-1] else 0
      # consume 1 if available
      if (stock[i] > 0) stock[i] <- stock[i] - 1
    }
    covered <- days[stock < lag(stock, default=0)]  # days where consumption occurred
    return(covered)
  }
}

# derive per-block features
derive_block_features <- function(cohort_ids, K=4, len=90, carry="cap30") {
  # precompute rx by class; labs; elig; death
  rx <- raw$rx %>% filter(pat_id %in% cohort_ids)
  labs <- raw$labs %>% filter(pat_id %in% cohort_ids)
  elig <- raw$elig %>% filter(pat_id %in% cohort_ids)
  death<- raw$death %>% filter(pat_id %in% cohort_ids)

  out <- map_dfr(cohort_ids, function(id){
    idx <- cohort %>% filter(pat_id==id) %>% pull(index_date)
    A0  <- cohort %>% filter(pat_id==id) %>% pull(A0)

    blks <- make_blocks(idx, K, len)
    # death / enrollment
    ddate <- death %>% filter(pat_id==id) %>% pull(death_date)
    e <- elig %>% filter(pat_id==id) %>% arrange(enr_start)

    # gather rx fills
    rxi <- rx %>% filter(pat_id==id)
    # split by class
    rxa <- rxi %>% filter(str_detect(drug_class,"ART_A"))
    rxb <- rxi %>% filter(str_detect(drug_class,"ART_B"))
    covA <- covered_days(rxa %>% select(fill_date, days_supply), carry = carry)
    covB <- covered_days(rxb %>% select(fill_date, days_supply), carry = carry)

    # labs
    li <- labs %>% filter(pat_id==id)

    # per-block measures
    map_dfr(1:K, function(k){
      s <- blks$start[k]; ebd <- blks$end[k]

      # enrollment alive through end
      alive <- ifelse(length(ddate)==0, TRUE, ebd < ddate)
      enrolled <- any(e$enr_start <= s & e$enr_end >= ebd)

      # coverage in block (A vs B)
      covA_k <- sum(covA >= s & covA <= ebd)
      covB_k <- sum(covB >= s & covB <= ebd)
      pdcA <- covA_k / 90; pdcB <- covB_k / 90
      pdc  <- pmin(1, pdcA + pdcB)  # total ART coverage (any)
      Ak   <- as.integer(pdcA >= pdcB) # "on A" majority coverage

      # monitoring: any VL in block
      M    <- any(li$test=="VL" & li$lab_date >= s & li$lab_date <= ebd)

      # simple L_k proxy: last observed VL category & utilization proxy (count labs/visits)
      last_vl <- li %>% filter(test=="VL", lab_date < s) %>% arrange(desc(lab_date)) %>% slice(1) %>% pull(result)
      Lk <- case_when(
        is.na(last_vl) ~ 0,
        last_vl < 200  ~ -1,
        last_vl < 1000 ~ 0.5,
        TRUE ~ 1.5
      )

      tibble(pat_id=id, k=k, start=s, end=ebd,
             A = Ak, PDC = pdc, L = Lk, M = as.integer(M),
             C = as.integer(enrolled & alive),
             D = as.integer(!alive & ebd >= ddate & s <= ddate))
    }) %>% mutate(A0 = A0)
  })
  out
}

K <- 4; LEN <- 90
blocks <- derive_block_features(cohort$pat_id, K=K, len=LEN, carry="cap30")
head(blocks)




# confirmed failure utility on raw labs
confirmed_failure <- function(labs_df, id, start, end) {
  li <- labs_df %>% filter(pat_id==id, test=="VL", lab_date >= start, lab_date <= end) %>%
    arrange(lab_date)
  if (nrow(li) < 2) return(0L)
  hi <- which(li$result >= 200)
  if (length(hi) < 2) return(0L)
  any(diff(li$lab_date[hi]) >= 7L) * 1L
}

# compute Y by horizon and per-block Y_k, D_{k-1}
mk_outcomes <- function(cohort_ids, K=4, len=90){
  labs <- raw$labs %>% filter(pat_id %in% cohort_ids)
  death<- raw$death %>% filter(pat_id %in% cohort_ids)

  res <- map_dfr(cohort_ids, function(id){
    idx <- cohort %>% filter(pat_id==id) %>% pull(index_date)
    blks <- make_blocks(idx, K, len)
    # Y (12m): confirmed failure by end of block K
    Y <- confirmed_failure(labs, id, blks$start[1], blks$end[K])

    # per-block cumulative incidence Y_k (LOCF)
    Yk <- integer(K); cum <- 0L
    for (k in 1:K) {
      if (cum==0L) {
        cum <- confirmed_failure(labs, id, blks$start[1], blks$end[k])
      }
      Yk[k] <- cum
    }
    # competing death indicator aligned as D_0..D_{K-1} (death before evaluating Y_{k+1})
    ddate <- death %>% filter(pat_id==id) %>% pull(death_date)
    Dcomp <- integer(K);  # will store D_0..D_{K-1} as length K (we'll drop last later)
    for (k in 0:(K-1)) {
      Dcomp[k+1] <- if (length(ddate)==0) 0L else as.integer(ddate <= blks$end[k+1])
    }
    tibble(pat_id=id, Y=Y,
           !!!set_names(as.list(Yk), paste0("Y_", 1:K)),
           !!!set_names(as.list(Dcomp[1:K]), paste0("D_", 0:(K-1))))
  })
  res
}

outcomes <- mk_outcomes(cohort$pat_id, K=K, len=LEN)
head(outcomes)


# wide pivot
wide <- blocks %>%
  select(pat_id, k, A, PDC, L, M, C, D, A0) %>%
  mutate(k = paste0("", k)) %>%
  pivot_wider(names_from = k,
              values_from = c(A,PDC,L,M,C,D),
              names_glue = "{.value}{k}") %>%
  left_join(outcomes, by="pat_id") %>%
  left_join(cohort %>% select(pat_id, sex, age), by="pat_id") %>%
  # sanity: bound PDC
  mutate(across(starts_with("PDC"), ~pmin(1, pmax(0, .x))))

# show columns
names(wide)


