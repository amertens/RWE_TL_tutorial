# Set seed for reproducibility
set.seed(12345)

# Load required libraries
library(survival)
library(concrete)

# Load and prepare the PBC dataset
# Select relevant columns and remove subjects with missing treatment info.
data <- survival::pbc[, c("time", "status", "trt", "age", "sex", "albumin")]
data <- subset(data, !is.na(data$trt))
# Recode treatment so that 0 indicates placebo and 1 indicates D-penicillamine
data$trt <- data$trt - 1

# Define target times for estimation (e.g., biannual time points between 3 and 6 years)
target_times <- 365.25/2 * (6:12)

# Format the analysis arguments using the concrete package's formatArguments() function
ConcreteArgs <- formatArguments(
  DataTable   = data,
  EventTime   = "time",
  EventType   = "status",
  Treatment   = "trt",
  Intervention= 0:1,              # Two static interventions: treat all with 0 and treat all with 1
  TargetTime  = target_times,     # Time points (in days) at which to estimate the risks
  TargetEvent = 1:2,              # Two competing events (e.g., 1 for transplant and 2 for death)
  CVArg       = list(V = 10),     # Use 10-fold cross-validation
  Verbose     = FALSE
)

# Run the one-step TMLE estimation using doConcrete()
ConcreteEst <- doConcrete(ConcreteArgs)

# Extract the output:
# Here, we request the risk difference (RD) estimand between the interventions,
# with simultaneous confidence intervals computed.
ConcreteOut <- getOutput(ConcreteEst, Estimand = "RD", Simultaneous = TRUE)

# Print the first few rows of the output table
print(head(ConcreteOut))

# Plot the estimated counterfactual risk curves and risk difference
# (The plot will include null reference lines for visual comparison.)
plot(ConcreteOut, NullLine = TRUE, ask = FALSE)
