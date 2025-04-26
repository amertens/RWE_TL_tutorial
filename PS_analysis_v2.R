############################################################
# 1) LOAD YOUR NEW DATA
############################################################

# Suppose you have "df_new" with these columns:
#   - treatment (0=non-SOF, 1=SOF)
#   - event (0/1 for AKI)
#   - observed_time (follow-up time)
#   - covariates: age, sex, race, region, cirrhosis, CKD, mental_ill, ...
#   - possibly an imputed version if you had missingness

df <- df_final   # <--- just rename to keep code consistent

############################################################
# 2) SPECIFY YOUR PROPENSITY SCORE FORMULA
############################################################

ps_formula <- as.formula(
  "treatment ~ age + sex + race + region + cirrhosis + hiv + mental_illness +
               ckd + diabetes + hypertension + baseline_gfr + bmi"
)
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
  include = c("observed_time", "event")
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
cox_matched <- coxph(Surv(time, event) ~ treatment, data=matched_df)
summary(cox_matched)

# Possibly, a fully adjusted Cox:
cox_adj <- coxph(Surv(time, event) ~ treatment + age + sex + race + region + cirrhosis + ckd + mental_illness + diabetes + hypertension + baseline_gfr + bmi,
                 data=matched_df)
summary(cox_adj)
