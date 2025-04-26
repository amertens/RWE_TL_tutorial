library(survtmle);  library(SuperLearner); set.seed(1)

n  <- 1000
W1 <- rnorm(n)                  # baseline covariate
A  <- rbinom(n, 1, 0.5)         # randomised treatment

## latent event & censoring time
lambda0 <- 0.08;  HR <- 0.6
Tlat <- rexp(n, lambda0 * ifelse(A==1, HR, 1))
Clat <- rexp(n, 0.02)

time  <- pmin(Tlat, Clat)
event <- as.integer(Tlat <= Clat)   # 1 = event, 0 = censored

ftime <- ceiling(time)          # <- make them integers
ftype <- event         # 0/1   (one type of failure)

sl_lib_f <- c("SL.glm","SL.mean")   # outcome hazard / iterated mean
sl_lib_c <- "SL.glm"                # censoring hazard
sl_lib_g <- "SL.glm"                # treatment model (not required in an RCT
# but keeps the template general)
t0  <- 5                    # years (can be a vector)

fit <- survtmle(
  ftime       = ftime,
  ftype       = ftype,      # 0 = censored, 1 = event
  trt         = A,
  adjustVars  = data.frame(W1 = W1),
  t0          = t0,

  ### nuisance-parameter learners
  SL.ftime    = sl_lib_f,   # failure model
  SL.ctime    = sl_lib_c,   # censoring model
  SL.trt      = sl_lib_g,   # treatment model (optional in RCT)

  method      = "hazard"    # canonical survival-TMLE (default)
)

fit$est                      # marginal survival estimates

S1 <- fit$est[1,1]
S0 <- fit$est[2,1]
Lambda1 <- -log(S1);         Lambda0 <- -log(S0)
cumHR   <- Lambda1 / Lambda0

## influence-curve–based SE & 95% CI for log-cHR
ic  <- (fit$ic[,2] - fit$ic[,1]) / (S1*log(S0/S1))
se  <- sd(ic)/sqrt(n)
ci  <- exp( log(cumHR) + c(-1,1)*1.96*se )

cat("Cumulative-hazard ratio =", round(cumHR,3),
    "; 95% CI:", paste(round(ci,3), collapse=" – "), "\n")
