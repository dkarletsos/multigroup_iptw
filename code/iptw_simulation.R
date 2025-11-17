################################################################################
# MUST RUN THIS BLOCK SEPARATELY TO SET CURRENT DIRECTORY
library(rstudioapi)



if (rstudioapi::isAvailable()) {
  # Get the directory of the currently active script
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  print(script_dir)
} else {
  stop("rstudioapi is not available. Make sure you're running this in RStudio.")
}
setwd(script_dir)
################################################################################

################################################################################
# Inverse Probability Weighting for Multiple Groups (Trial Membership)
# Simulation Study for IPD Meta-Analysis
# R Implementation
################################################################################

# Clear workspace
rm(list = ls())

# Set working directory (modify as needed)
# setwd("path/to/your/directory")

################################################################################
# INSTALL AND LOAD REQUIRED PACKAGES
################################################################################

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Required packages
required_packages <- c("nnet", "survival", "ggplot2", "dplyr", "tidyr", 
                       "gridExtra", "scales", "tableone")

install_if_missing(required_packages)

# Set seed for reproducibility
set.seed(12345)

cat("\n")
cat("================================================================================\n")
cat("IPTW FOR MULTIPLE GROUPS SIMULATION - R IMPLEMENTATION\n")
cat("================================================================================\n\n")

################################################################################
# PART 1: DATA GENERATION FUNCTION
################################################################################

cat("================================================================================\n")
cat("PART 1: GENERATING IPD DATA WITH HETEROGENEOUS TRIALS\n")
cat("================================================================================\n\n")

generate_ipd_data <- function(n_trials = 5, n_per_trial = 300) {
  #' Generate IPD with heterogeneous baseline characteristics across trials
  #' 
  #' @param n_trials Number of trials to simulate
  #' @param n_per_trial Number of participants per trial
  #' @return Data frame with IPD from all trials
  
  # Trial-specific parameters (simulating different populations)
  trial_params <- list(
    # Trial 1: Younger, lower severity
    list(age_mean = 55, age_sd = 10, severity_mean = 45, 
         severity_sd = 15, female_prop = 0.45),
    # Trial 2: Older, higher severity
    list(age_mean = 70, age_sd = 8, severity_mean = 68, 
         severity_sd = 12, female_prop = 0.55),
    # Trial 3: Moderate characteristics
    list(age_mean = 62, age_sd = 12, severity_mean = 55, 
         severity_sd = 18, female_prop = 0.50),
    # Trial 4: Younger with higher severity
    list(age_mean = 58, age_sd = 11, severity_mean = 62, 
         severity_sd = 14, female_prop = 0.48),
    # Trial 5: Older with lower severity
    list(age_mean = 67, age_sd = 9, severity_mean = 48, 
         severity_sd = 16, female_prop = 0.52)
  )
  
  # Initialize empty list to store trial data
  trial_list <- list()
  
  # Generate data for each trial
  for (i in 1:n_trials) {
    params <- trial_params[[i]]
    
    # Generate correlated baseline covariates
    age <- rnorm(n_per_trial, mean = params$age_mean, sd = params$age_sd)
    
    # Disease severity correlated with age
    severity <- params$severity_mean + 0.3 * (age - params$age_mean) + 
                rnorm(n_per_trial, 0, params$severity_sd)
    severity <- pmax(0, pmin(100, severity))  # Bound between 0-100
    
    # Sex (female = 1)
    female <- rbinom(n_per_trial, 1, params$female_prop)
    
    # Comorbidity count (0-5) - related to age
    age_scaled <- (age - mean(age)) / sd(age)
    comorbid_prob <- plogis(-0.5 + 0.5 * age_scaled + rnorm(n_per_trial, 0, 0.3))
    comorbidities <- rbinom(n_per_trial, 5, comorbid_prob)
    
    # Generate survival outcomes
    # log(hazard) depends on baseline covariates
    linear_predictor <- -3.5 + 
                       0.03 * (age - 60) +           # Age effect (centered)
                       0.02 * (severity - 50) +      # Severity effect (centered)
                       0.15 * comorbidities +        # Comorbidity effect
                       -0.2 * female                 # Female protective effect
    
    # Generate survival times from Weibull distribution
    shape <- 1.2  # Weibull shape parameter
    scale <- exp(-linear_predictor / shape)
    u <- runif(n_per_trial)
    time_to_event <- scale * (-log(u))^(1/shape)
    
    # Administrative censoring at 5 years
    censor_time <- 5
    time <- pmin(time_to_event, censor_time)
    event <- as.numeric(time_to_event <= censor_time)
    
    # Create trial dataset (PLACEBO ARM ONLY for this analysis)
    trial_data <- data.frame(
      trial_id = as.factor(i),
      patient_id = paste0("T", i, "_", sprintf("%03d", 1:n_per_trial)),
      age = age,
      female = female,
      severity = severity,
      comorbidities = comorbidities,
      time = time,
      event = event,
      stringsAsFactors = FALSE
    )
    
    trial_list[[i]] <- trial_data
  }
  
  # Combine all trials
  ipd_data <- bind_rows(trial_list)
  
  return(ipd_data)
}

# Generate the IPD dataset
ipd_data <- generate_ipd_data(n_trials = 5, n_per_trial = 300)

cat("IPD Data Generated\n")
cat("Total sample size:", nrow(ipd_data), "\n")
cat("Number of trials:", length(unique(ipd_data$trial_id)), "\n\n")

cat("Sample size per trial:\n")
print(table(ipd_data$trial_id))
cat("\n")

################################################################################
# PART 2: ASSESS PRE-WEIGHTING BALANCE
################################################################################

cat("================================================================================\n")
cat("PART 2: PRE-WEIGHTING BALANCE ASSESSMENT\n")
cat("================================================================================\n\n")

# Function to calculate standardized mean difference (SMD)
calculate_smd <- function(data, var, group_var) {
  #' Calculate SMD comparing each group to overall mean
  #' 
  #' @param data Data frame
  #' @param var Variable name to calculate SMD for
  #' @param group_var Grouping variable name
  #' @return Data frame with SMDs by group
  
  overall_mean <- mean(data[[var]], na.rm = TRUE)
  overall_sd <- sd(data[[var]], na.rm = TRUE)
  
  smds <- data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      group_mean = mean(.data[[var]], na.rm = TRUE),
      smd = (group_mean - overall_mean) / overall_sd,
      .groups = 'drop'
    )
  
  return(smds)
}

# Function to calculate weighted SMD
calculate_weighted_smd <- function(data, var, group_var, weight_var) {
  #' Calculate weighted SMD comparing each group to weighted overall mean
  #' 
  #' @param data Data frame
  #' @param var Variable name to calculate SMD for
  #' @param group_var Grouping variable name
  #' @param weight_var Weight variable name
  #' @return Data frame with weighted SMDs by group
  
  # Weighted overall mean
  overall_mean <- weighted.mean(data[[var]], data[[weight_var]], na.rm = TRUE)
  
  # Weighted overall SD
  weighted_var <- sum(data[[weight_var]] * (data[[var]] - overall_mean)^2, 
                      na.rm = TRUE) / sum(data[[weight_var]])
  overall_sd <- sqrt(weighted_var)
  
  smds <- data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      weighted_mean = weighted.mean(.data[[var]], .data[[weight_var]], na.rm = TRUE),
      smd = (weighted_mean - overall_mean) / overall_sd,
      .groups = 'drop'
    )
  
  return(smds)
}

# Calculate pre-weighting SMDs
vars_to_check <- c("age", "female", "severity", "comorbidities")

cat("PRE-WEIGHTING STANDARDIZED MEAN DIFFERENCES:\n\n")

pre_weight_smds <- list()
for (var in vars_to_check) {
  smd_results <- calculate_smd(ipd_data, var, "trial_id")
  pre_weight_smds[[var]] <- smd_results
  
  cat(paste0("Variable: ", toupper(var), "\n"))
  print(smd_results, n = Inf)
  cat(paste0("Max absolute SMD: ", round(max(abs(smd_results$smd)), 3), "\n\n"))
}

# Summary statistics by trial
cat("BASELINE CHARACTERISTICS BY TRIAL (UNWEIGHTED):\n\n")

baseline_summary <- ipd_data %>%
  group_by(trial_id) %>%
  summarise(
    N = n(),
    Age_mean = round(mean(age), 1),
    Age_sd = round(sd(age), 1),
    Female_pct = round(100 * mean(female), 1),
    Severity_mean = round(mean(severity), 1),
    Severity_sd = round(sd(severity), 1),
    Comorbid_mean = round(mean(comorbidities), 2),
    Comorbid_sd = round(sd(comorbidities), 2),
    Events = sum(event),
    Event_rate = round(100 * mean(event), 1),
    .groups = 'drop'
  )

print(baseline_summary, n = Inf)
cat("\n")

################################################################################
# PART 3: ESTIMATE PROPENSITY SCORES (TRIAL MEMBERSHIP)
################################################################################

cat("================================================================================\n")
cat("PART 3: ESTIMATING TRIAL MEMBERSHIP PROPENSITY SCORES\n")
cat("================================================================================\n\n")

# Fit multinomial logistic regression
# Outcome: trial membership (trial_id)
# Predictors: baseline covariates
cat("Fitting multinomial logistic regression model...\n")

ps_model <- multinom(
  trial_id ~ age + female + severity + comorbidities, 
  data = ipd_data, 
  trace = FALSE
)

cat("Multinomial logistic regression model fitted.\n\n")

cat("Model summary:\n")
print(summary(ps_model))
cat("\n")

# Predict propensity scores (probability of being in each trial)
ps_matrix <- predict(ps_model, type = "probs")

# Ensure ps_matrix is a matrix (it's a vector if only 2 trials)
if (!is.matrix(ps_matrix)) {
  ps_matrix <- cbind(ps_matrix, 1 - ps_matrix)
}

# For each participant, extract their probability of being in their observed trial
ipd_data$ps <- sapply(1:nrow(ipd_data), function(i) {
  trial <- as.numeric(ipd_data$trial_id[i])
  ps_matrix[i, trial]
})

cat("PROPENSITY SCORE DIAGNOSTICS:\n\n")
cat("Propensity score summary:\n")
print(summary(ipd_data$ps))
cat("\n")

cat("Propensity scores by trial:\n")
ps_by_trial <- ipd_data %>% 
  group_by(trial_id) %>% 
  summarise(
    Mean_PS = round(mean(ps), 3),
    Median_PS = round(median(ps), 3),
    Min_PS = round(min(ps), 3),
    Max_PS = round(max(ps), 3),
    .groups = 'drop'
  )
print(ps_by_trial, n = Inf)
cat("\n")

################################################################################
# PART 4: CALCULATE STABILIZED IPTW WEIGHTS
################################################################################

cat("================================================================================\n")
cat("PART 4: CALCULATING STABILIZED IPTW WEIGHTS\n")
cat("================================================================================\n\n")

# Calculate marginal probability for each trial (proportion in each trial)
marginal_probs <- ipd_data %>%
  group_by(trial_id) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(marginal_prob = n / sum(n))

# Add marginal probabilities to data
ipd_data <- ipd_data %>%
  left_join(marginal_probs %>% select(trial_id, marginal_prob), 
            by = "trial_id")

# Calculate stabilized weights: sw = P(Trial = j) / P(Trial = j | X)
ipd_data$sw <- ipd_data$marginal_prob / ipd_data$ps

# Calculate unstabilized weights for comparison: w = 1 / P(Trial = j | X)
ipd_data$uw <- 1 / ipd_data$ps

cat("Stabilized weights calculated.\n\n")

cat("WEIGHT DIAGNOSTICS:\n\n")
cat("Stabilized weights summary:\n")
print(summary(ipd_data$sw))
cat("\n")

cat("Unstabilized weights summary:\n")
print(summary(ipd_data$uw))
cat("\n")

cat("Stabilized weights by trial:\n")
sw_by_trial <- ipd_data %>% 
  group_by(trial_id) %>% 
  summarise(
    Mean_SW = round(mean(sw), 3),
    Median_SW = round(median(sw), 3),
    Min_SW = round(min(sw), 3),
    Max_SW = round(max(sw), 3),
    SD_SW = round(sd(sw), 3),
    .groups = 'drop'
  )
print(sw_by_trial, n = Inf)
cat("\n")

# Calculate effective sample size
ess_stabilized <- sum(ipd_data$sw)^2 / sum(ipd_data$sw^2)
ess_unstabilized <- sum(ipd_data$uw)^2 / sum(ipd_data$uw^2)

cat("EFFECTIVE SAMPLE SIZE:\n")
cat("Original sample size:", nrow(ipd_data), "\n")
cat("Effective sample size (stabilized weights):", round(ess_stabilized, 0), 
    "(", round(100 * ess_stabilized / nrow(ipd_data), 1), "%)\n")
cat("Effective sample size (unstabilized weights):", round(ess_unstabilized, 0), 
    "(", round(100 * ess_unstabilized / nrow(ipd_data), 1), "%)\n\n")

################################################################################
# PART 5: ASSESS POST-WEIGHTING BALANCE
################################################################################

cat("================================================================================\n")
cat("PART 5: POST-WEIGHTING BALANCE ASSESSMENT\n")
cat("================================================================================\n\n")

cat("POST-WEIGHTING STANDARDIZED MEAN DIFFERENCES:\n\n")

post_weight_smds <- list()
for (var in vars_to_check) {
  smd_results <- calculate_weighted_smd(ipd_data, var, "trial_id", "sw")
  post_weight_smds[[var]] <- smd_results
  
  cat(paste0("Variable: ", toupper(var), "\n"))
  print(smd_results, n = Inf)
  cat(paste0("Max absolute SMD: ", round(max(abs(smd_results$smd)), 3), "\n\n"))
}

# Create comparison table
cat("SMD COMPARISON: PRE vs POST WEIGHTING:\n\n")

smd_comparison <- data.frame()
for (var in vars_to_check) {
  pre <- pre_weight_smds[[var]]
  post <- post_weight_smds[[var]]
  
  for (i in 1:nrow(pre)) {
    improvement <- ifelse(abs(post$smd[i]) < abs(pre$smd[i]), "Improved", "No change")
    smd_comparison <- rbind(smd_comparison, data.frame(
      Variable = var,
      Trial = as.character(pre$trial_id[i]),
      SMD_Unweighted = round(pre$smd[i], 3),
      SMD_Weighted = round(post$smd[i], 3),
      Improvement = improvement,
      stringsAsFactors = FALSE
    ))
  }
}

print(smd_comparison, row.names = FALSE)
cat("\n")

# Summary of balance improvement
cat("BALANCE IMPROVEMENT SUMMARY:\n")
n_pre_imbalanced <- sum(abs(smd_comparison$SMD_Unweighted) > 0.1)
n_post_imbalanced <- sum(abs(smd_comparison$SMD_Weighted) > 0.1)
pct_improved <- 100 * mean(smd_comparison$Improvement == "Improved")

cat("Number of comparisons with |SMD| > 0.1 before weighting:", n_pre_imbalanced, "\n")
cat("Number of comparisons with |SMD| > 0.1 after weighting:", n_post_imbalanced, "\n")
cat("Proportion of comparisons improved:", round(pct_improved, 1), "%\n\n")

################################################################################
# PART 6: VISUALIZATIONS
################################################################################

cat("================================================================================\n")
cat("PART 6: CREATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Prepare data for love plot
love_plot_data <- smd_comparison %>%
  pivot_longer(cols = c(SMD_Unweighted, SMD_Weighted),
               names_to = "Weighting",
               values_to = "SMD") %>%
  mutate(
    Weighting = factor(Weighting, 
                      levels = c("SMD_Unweighted", "SMD_Weighted"),
                      labels = c("Unweighted", "Weighted")),
    Variable_Trial = paste0(Variable, " (Trial ", Trial, ")")
  )

# Create love plot
p1 <- ggplot(love_plot_data, aes(x = abs(SMD), y = Variable_Trial, 
                                  color = Weighting, shape = Weighting)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", alpha = 0.6) +
  scale_color_manual(values = c("Unweighted" = "darkblue", "Weighted" = "darkgreen")) +
  scale_shape_manual(values = c("Unweighted" = 16, "Weighted" = 17)) +
  labs(
    title = "Covariate Balance Before and After IPTW",
    subtitle = "Standardized Mean Differences by Trial",
    x = "Absolute Standardized Mean Difference",
    y = "",
    color = "",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

ggsave("love_plot.png", p1, width = 10, height = 8, dpi = 300)
cat("Love plot saved: love_plot.png\n")

# Weight distribution plot
p2 <- ggplot(ipd_data, aes(x = sw)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Stabilized IPTW Weights",
    x = "Stabilized Weight",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

ggsave("weight_distribution.png", p2, width = 8, height = 6, dpi = 300)
cat("Weight distribution plot saved: weight_distribution.png\n")

# Propensity score distribution by trial
p3 <- ggplot(ipd_data, aes(x = ps, fill = trial_id)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Propensity Score Distribution by Trial",
    x = "Propensity Score (Probability of Observed Trial Membership)",
    y = "Density",
    fill = "Trial"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

ggsave("propensity_score_distribution.png", p3, width = 10, height = 6, dpi = 300)
cat("Propensity score distribution plot saved: propensity_score_distribution.png\n\n")

################################################################################
# PART 7: SURVIVAL ANALYSIS (DEMONSTRATION)
################################################################################

cat("================================================================================\n")
cat("PART 7: SURVIVAL ANALYSIS DEMONSTRATION\n")
cat("================================================================================\n\n")

# Add treatment arm for demonstration (50% treatment, 50% placebo)
set.seed(123)
ipd_data$treatment <- as.numeric(as.numeric(rownames(ipd_data)) %% 2)

# Modify survival times for treatment arm (true HR = 0.70)
treatment_effect <- log(0.70)
ipd_data$time_modified <- ifelse(
  ipd_data$treatment == 1,
  ipd_data$time * exp(-treatment_effect),
  ipd_data$time
)
ipd_data$event_modified <- ifelse(
  ipd_data$treatment == 1 & ipd_data$time_modified > 5,
  0, 
  ipd_data$event
)
ipd_data$time_modified <- pmin(ipd_data$time_modified, 5)

# Unweighted Cox model
cat("UNWEIGHTED ANALYSIS:\n")
cox_unweighted <- coxph(Surv(time_modified, event_modified) ~ treatment, 
                        data = ipd_data)
print(summary(cox_unweighted))
cat("\n")

# Weighted Cox model (using stabilized weights)
cat("WEIGHTED ANALYSIS (Stabilized IPTW):\n")
cox_weighted <- coxph(Surv(time_modified, event_modified) ~ treatment, 
                      data = ipd_data,
                      weights = sw,
                      robust = TRUE)
print(summary(cox_weighted))
cat("\n")

# Extract results
hr_unweighted <- exp(coef(cox_unweighted))
ci_unweighted <- exp(confint(cox_unweighted))
hr_weighted <- exp(coef(cox_weighted))
ci_weighted <- exp(confint(cox_weighted))

cat("COMPARISON OF TREATMENT EFFECT ESTIMATES:\n")
cat("Unweighted HR:", round(hr_unweighted, 3), 
    "95% CI: (", round(ci_unweighted[1], 3), "-", round(ci_unweighted[2], 3), ")\n")
cat("Weighted HR:", round(hr_weighted, 3), 
    "95% CI: (", round(ci_weighted[1], 3), "-", round(ci_weighted[2], 3), ")\n\n")

################################################################################
# PART 8: SAVE RESULTS
################################################################################

cat("================================================================================\n")
cat("PART 8: SAVING RESULTS\n")
cat("================================================================================\n\n")

# Save processed dataset
write.csv(ipd_data, "ipd_data_with_weights.csv", row.names = FALSE)
cat("IPD data with weights saved: ipd_data_with_weights.csv\n")

# Save SMD comparison table
write.csv(smd_comparison, "smd_comparison.csv", row.names = FALSE)
cat("SMD comparison table saved: smd_comparison.csv\n")

# Save baseline characteristics table
write.csv(baseline_summary, "baseline_characteristics.csv", row.names = FALSE)
cat("Baseline characteristics table saved: baseline_characteristics.csv\n")

# Create summary report
summary_report <- data.frame(
  Metric = c(
    "Original Sample Size",
    "Effective Sample Size (Stabilized)",
    "ESS Retention (%)",
    "Max Pre-Weight |SMD|",
    "Max Post-Weight |SMD|",
    "Comparisons |SMD| > 0.1 (Pre)",
    "Comparisons |SMD| > 0.1 (Post)",
    "Unweighted HR",
    "Weighted HR"
  ),
  Value = c(
    nrow(ipd_data),
    round(ess_stabilized, 0),
    round(100 * ess_stabilized / nrow(ipd_data), 1),
    round(max(abs(smd_comparison$SMD_Unweighted)), 3),
    round(max(abs(smd_comparison$SMD_Weighted)), 3),
    n_pre_imbalanced,
    n_post_imbalanced,
    round(hr_unweighted, 3),
    round(hr_weighted, 3)
  )
)

write.csv(summary_report, "summary_report.csv", row.names = FALSE)
cat("Summary report saved: summary_report.csv\n\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("================================================================================\n")
cat("SIMULATION COMPLETE!\n")
cat("================================================================================\n\n")

cat("All results have been saved:\n")
cat("  - ipd_data_with_weights.csv\n")
cat("  - smd_comparison.csv\n")
cat("  - baseline_characteristics.csv\n")
cat("  - summary_report.csv\n")
cat("  - love_plot.png\n")
cat("  - weight_distribution.png\n")
cat("  - propensity_score_distribution.png\n\n")

cat("Key findings:\n")
cat("  - Sample size: N =", nrow(ipd_data), "\n")
cat("  - Effective sample size:", round(ess_stabilized, 0), 
    "(", round(100 * ess_stabilized / nrow(ipd_data), 1), "%)\n")
cat("  - Max pre-weighting |SMD|:", round(max(abs(smd_comparison$SMD_Unweighted)), 3), "\n")
cat("  - Max post-weighting |SMD|:", round(max(abs(smd_comparison$SMD_Weighted)), 3), "\n")
cat("  - Imbalanced comparisons reduced from", n_pre_imbalanced, "to", n_post_imbalanced, "\n\n")

cat("================================================================================\n")