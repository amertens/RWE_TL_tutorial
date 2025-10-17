############################################################
# 1) LOAD YOUR NEW DATA
############################################################

# Suppose you have "df_new" with these columns:
#   - treatment (0=non-SOF, 1=SOF)
#   - event (0/1 for AKI)
#   - observed_time (follow-up time)
#   - covariates: age, sex, race, region, cirrhosis, CKD, mental_ill, ...
#   - possibly an imputed version if you had missingness


df <- read_csv(here::here("data/sim_hcv_aki_complex.csv"))   # <--- just rename to keep code consistent
head(df)

df$id
colnames(df)
df$follow_time

############################################################
# 2) SPECIFY YOUR PROPENSITY SCORE FORMULA
############################################################

ps_formula <- as.formula(
  "treatment ~ age + sex_male + race + region + cirrhosis + hiv +
               ckd + diabetes + hypertension + sepsis + bmi"
)

# [8] "heart_failure"    "sepsis"           "dehydration"      "obstruction"      "cirrhosis"        "portal_htn"       "esld"
# [15] "hiv"              "diabetes"         "hypertension"     "bmi"              "overweight_obese" "smoking"          "alcohol"
# [22] "substance_abuse"  "cancer"           "chemo"            "nsaid"            "acearb"           "diuretic"         "aminoglycoside"
# [29] "contrast"         "statin"           "aspirin"          "beta_blocker"     "ccb"              "art"              "treatment"
# [36] "switch"           "follow_time"      "event"

# Example logistic regression for propensity:
ps_model <- glm(ps_formula, data = df, family = binomial())

df$ps <- predict(ps_model, type = "response")

############################################################
# 3) MATCH ON PROPENSITY SCORE
############################################################
library(MatchIt)


matchit_out <- matchit(
  ps_formula,
  data = df,
  method = "nearest",
  caliper = 0.2,
  ratio = 1
)

matched_df <- match.data(
  matchit_out,
  data = df,  # re-specify the full data
  include = c("follow_time", "event")
)


############################################################
# 4) CHECK BALANCE, ASSESS OUTCOME
############################################################

# Covariate balance
library(cobalt)
bal.tab(matchit_out, un = TRUE, m.threshold = 0.1)

# Summaries
library(tableone)
vars_to_check <- c("age","sex","race","region","cirrhosis","ckd","mental_illness",
                   "diabetes","hypertension","baseline_gfr","bmi")
table_matched <- CreateTableOne(data=matched_df, vars=vars_to_check,
                                strata="treatment")
print(table_matched, smd=TRUE)

# Estimate effect
# For instance, a Cox model on matched sample
library(survival)
cox_matched <- coxph(Surv(follow_time, event) ~ treatment, data=matched_df)
summary(cox_matched)

# Possibly, a fully adjusted Cox:
cox_adj <- coxph(Surv(follow_time, event) ~ treatment + age + sex_male + race + region + cirrhosis + ckd  + diabetes + hypertension  + bmi,
                 data=matched_df)
summary(cox_adj)
