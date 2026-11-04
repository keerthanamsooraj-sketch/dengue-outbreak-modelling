#INSTALL AND LOAD PACKAGES
packages <- c("outbreaks", "deSolve")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}
library(outbreaks)
library(deSolve)
#IMPORT ACTUAL DENGUE OUTBREAK DATA
data("dengue_yap_2011", package = "outbreaks")
dengue <- dengue_yap_2011
# Inspect dataset
head(dengue)
str(dengue)
summary(dengue)
# Check missing values
colSums(is.na(dengue))

#PREPARE DAILY CASE DATA
# Convert days into a continuous daily sequence
dengue$day <- as.integer(
  dengue$nr - min(dengue$nr, na.rm = TRUE)
)
# Combine any duplicate days
daily_cases <- aggregate(
  value ~ day,
  data = dengue,
  FUN = sum
)
# Include days with zero reported cases
all_days <- data.frame(
  day = 0:max(daily_cases$day)
)

daily <- merge(
  all_days,
  daily_cases,
  by = "day",
  all.x = TRUE
)

daily$value[is.na(daily$value)] <- 0

# Assign dates
start_date <- min(dengue$onset_date, na.rm = TRUE)

daily$date <- start_date + daily$day

# Calculate cumulative cases
daily$cumulative_cases <- cumsum(daily$value)

# Check data
head(daily)
tail(daily)

cat("Total reported cases:", sum(daily$value), "\n")

#DESCRIPTIVE EPIDEMIOLOGY
population <- 7391

total_cases <- sum(daily$value)

attack_rate <- (total_cases / population) * 100

peak_cases <- max(daily$value)

peak_date <- daily$date[
  which.max(daily$value)
]

cat("Total cases:", total_cases, "\n")
cat("Attack rate (%):", round(attack_rate, 2), "\n")
cat("Maximum daily cases:", peak_cases, "\n")
cat("Date of maximum daily cases:", as.character(peak_date), "\n")

# VISUALISE THE EPIDEMIC CURVE
plot(
  daily$date,
  daily$value,
  type = "h",
  lwd = 2,
  col = "steelblue",
  xlab = "Date",
  ylab = "Daily reported cases",
  main = "Dengue Epidemic Curve: Yap Islands, 2011-2012"
)

#CUMULATIVE CASE CURVE
plot(
  daily$date,
  daily$cumulative_cases,
  type = "l",
  lwd = 2,
  col = "darkred",
  xlab = "Date",
  ylab = "Cumulative reported cases",
  main = "Cumulative Dengue Cases"
)

#DEFINE SIR MODEL
sir_model <- function(time, state, parameters) {

  with(as.list(c(state, parameters)), {

    # New infection flow
    infection <- beta * S * I / N

    # Recovery/removal flow
    recovery <- gamma * I

    # Differential equations
    dS <- -infection

    dI <- infection - recovery

    dR <- recovery

    dC <- infection

    list(c(dS, dI, dR, dC))
  })
}

#PREPARE MODEL INPUTS
N <- 7391

gamma <- 1 / 5

times <- 0:max(daily$day)

observed_cumulative <- daily$cumulative_cases

# 9. DEFINE PARAMETER ESTIMATION FUNCTION
fit_sir <- function(gamma_value) {

  objective <- function(par) {

    # Estimate parameters
    beta <- exp(par[1])

    I0 <- exp(par[2])

    reporting_scale <- plogis(par[3])

    # Prevent invalid initial conditions
    if (I0 >= N) return(1e12)

    initial_state <- c(
      S = N - I0,
      I = I0,
      R = 0,
      C = 0
    )

    parameters <- c(
      beta = beta,
      gamma = gamma_value,
      N = N
    )

    # Solve differential equations
    result <- try(
      deSolve::ode(
        y = initial_state,
        times = times,
        func = sir_model,
        parms = parameters
      ),
      silent = TRUE
    )

    if (inherits(result, "try-error")) {
      return(1e12)
    }

    # Predicted cumulative reported cases
    predicted_cumulative <- reporting_scale *
      (I0 + result[, "C"])

    # Compare observed and predicted curves
    error <- sum(
      (log1p(observed_cumulative) -
       log1p(predicted_cumulative))^2
    )

    return(error)
  }

  # Optimise model parameters
  optim(
    par = c(log(0.3), log(1), qlogis(0.5)),
    fn = objective,
    method = "L-BFGS-B",
    lower = c(log(0.001), log(0.1), qlogis(0.001)),
    upper = c(log(5), log(1000), qlogis(0.999))
  )
}


# FIT MODEL


model_fit <- fit_sir(gamma)

# Extract estimated parameters
beta_est <- exp(model_fit$par[1])

initial_infected <- exp(model_fit$par[2])

reporting_scale <- plogis(model_fit$par[3])

R0_effective <- beta_est / gamma

cat("Estimated beta:", beta_est, "\n")

cat("Recovery rate:", gamma, "\n")

cat("Estimated initial infected:", initial_infected, "\n")

cat("Estimated reporting scale:", reporting_scale, "\n")

cat("Effective SIR reproduction metric:", R0_effective, "\n")

cat("Optimisation convergence code:", model_fit$convergence, "\n")


# GENERATE MODEL PREDICTIONS
initial_state <- c(
  S = N - initial_infected,
  I = initial_infected,
  R = 0,
  C = 0
)

parameters <- c(
  beta = beta_est,
  gamma = gamma,
  N = N
)

model_output <- as.data.frame(
  ode(
    y = initial_state,
    times = times,
    func = sir_model,
    parms = parameters
  )
)

# Calculate predicted cumulative reported cases
predicted_cumulative <- reporting_scale *
  (initial_infected + model_output$C)

# Calculate predicted daily reported cases
predicted_daily <- c(
  predicted_cumulative[1],
  diff(predicted_cumulative)
)

# Combine observed and predicted values
comparison <- data.frame(
  date = daily$date,
  observed = daily$value,
  predicted = predicted_daily,
  observed_cumulative = daily$cumulative_cases,
  predicted_cumulative = predicted_cumulative
)

head(comparison)

# PLOT OBSERVED VS PREDICTED DAILY CASES
plot(
  comparison$date,
  comparison$observed,
  type = "l",
  col = "black",
  lwd = 2,
  xlab = "Date",
  ylab = "Daily cases",
  main = "Observed vs Predicted Dengue Cases"
)

lines(
  comparison$date,
  comparison$predicted,
  col = "red",
  lwd = 2
)

legend(
  "topright",
  legend = c("Observed cases", "SIR predicted cases"),
  col = c("black", "red"),
  lty = 1,
  lwd = 2
)


# PLOT CUMULATIVE MODEL FIT
plot(
  comparison$date,
  comparison$observed_cumulative,
  type = "l",
  col = "black",
  lwd = 2,
  xlab = "Date",
  ylab = "Cumulative cases",
  main = "Observed vs Predicted Cumulative Cases"
)

lines(
  comparison$date,
  comparison$predicted_cumulative,
  col = "blue",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Observed", "Predicted"),
  col = c("black", "blue"),
  lty = 1,
  lwd = 2
)


# MODEL PERFORMANCE
# Root Mean Square Error
RMSE <- sqrt(
  mean(
    (comparison$observed_cumulative -
     comparison$predicted_cumulative)^2
  )
)

# Mean Absolute Error
MAE <- mean(
  abs(
    comparison$observed_cumulative -
    comparison$predicted_cumulative
  )
)

# Mean Absolute Percentage Error
# Exclude zero observed values
nonzero <- comparison$observed_cumulative > 0

MAPE <- mean(
  abs(
    (comparison$observed_cumulative[nonzero] -
     comparison$predicted_cumulative[nonzero]) /
      comparison$observed_cumulative[nonzero]
  )
) * 100

cat("RMSE:", RMSE, "\n")
cat("MAE:", MAE, "\n")
cat("MAPE:", MAPE, "%\n")

# SENSITIVITY ANALYSIS
infectious_periods <- c(3, 5, 7)

sensitivity_results <- do.call(
  rbind,
  lapply(infectious_periods, function(period) {

    gamma_value <- 1 / period

    fit <- fit_sir(gamma_value)

    beta_value <- exp(fit$par[1])

    data.frame(
      infectious_period = period,
      gamma = gamma_value,
      beta = beta_value,
      effective_R0 = beta_value / gamma_value,
      convergence = fit$convergence
    )
  })
)

print(sensitivity_results)



# EXPORT RESULTS
dir.create("Dengue_Portfolio", showWarnings = FALSE)

write.csv(
  comparison,
  "Dengue_Portfolio/dengue_model_comparison.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_results,
  "Dengue_Portfolio/sensitivity_analysis.csv",
  row.names = FALSE
)

# Save epidemic curve
png(
  "Dengue_Portfolio/epidemic_curve.png",
  width = 1200,
  height = 800
)

plot(
  daily$date,
  daily$value,
  type = "h",
  col = "steelblue",
  lwd = 2,
  xlab = "Date",
  ylab = "Daily cases",
  main = "Dengue Epidemic Curve: Yap Islands"
)

dev.off()

# Save model comparison graph
png(
  "Dengue_Portfolio/model_fit.png",
  width = 1200,
  height = 800
)

plot(
  comparison$date,
  comparison$observed_cumulative,
  type = "l",
  col = "black",
  lwd = 2,
  xlab = "Date",
  ylab = "Cumulative cases",
  main = "Observed vs Predicted Dengue Cases"
)

lines(
  comparison$date,
  comparison$predicted_cumulative,
  col = "blue",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Observed", "Predicted"),
  col = c("black", "blue"),
  lty = 1,
  lwd = 2
)

dev.off()

cat("Portfolio results exported successfully.\n")
cat("Population:", 7391, "\n")
cat("Number of observed days:", nrow(daily), "\n")
