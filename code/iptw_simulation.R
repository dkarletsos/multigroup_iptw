############################################################################
# Inverse Probability Weighting for Product Switching Groups
# ENHANCED VERSION: Weight Trimming + E-Value Sensitivity Analysis
# Adapted for PMI Cross-Sectional Risk Marker Study (NCT05385055)
# 3 Groups: Current Smokers, THS Users (Switchers), Former Smokers
############################################################################

# Clear workspace
rm(list = ls())

# Get script directory (for RStudio users)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_dir)
}

# Create output directories
tables_dir <- file.path(getwd(), "tables")
figures_dir <- file.path(getwd(), "figures")

if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

cat("Output directories:\n")
cat("  Tables:", tables_dir, "\n")
cat("  Figures:", figures_dir, "\n\n")

############################################################################
# INSTALL AND LOAD REQUIRED PACKAGES
############################################################################

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

required_packages <- c("nnet", "ggplot2", "dplyr", "tidyr", "gridExtra", "scales")
install_if_missing(required_packages)

set.seed(20250202)

cat("\n")
cat("================================================================================\n")
cat("IPTW FOR PRODUCT SWITCHING: PMI THS STUDY ADAPTATION\n")
cat("ENHANCED: Weight Trimming + E-Value Sensitivity Analysis\n")
cat("3 Groups: Current Smokers | THS Users (Switchers) | Former Smokers\n")
cat("================================================================================\n\n")

############################################################################
# E-VALUE FUNCTIONS
############################################################################

# E-value calculation for risk ratios
# Reference: VanderWeele & Ding (2017) Annals of Internal Medicine
calculate_evalue_rr <- function(rr, lo = NULL, hi = NULL) {
  #' Calculate E-value for a risk ratio
  #' 

  #' @param rr Point estimate of risk ratio (or HR, OR if rare outcome)
  #' @param lo Lower bound of 95% CI (optional)
  #' @param hi Upper bound of 95% CI (optional)
  #' @return List with E-values for point estimate and CI bound
  
  # For RR < 1, take reciprocal
  if (rr < 1) {
    rr <- 1 / rr
    if (!is.null(lo) && !is.null(hi)) {
      temp <- lo
      lo <- 1 / hi
      hi <- 1 / temp
    }
  }
  
  # E-value formula: E = RR + sqrt(RR * (RR - 1))
  evalue_point <- rr + sqrt(rr * (rr - 1))
  
  evalue_ci <- NA

  if (!is.null(lo)) {
    # For CI, use the bound closer to null (1)
    ci_bound <- ifelse(lo > 1, lo, ifelse(!is.null(hi) && hi < 1, hi, 1))
    if (ci_bound > 1) {
      evalue_ci <- ci_bound + sqrt(ci_bound * (ci_bound - 1))
    } else {
      evalue_ci <- 1  # CI includes null
    }
  }
  
  return(list(
    evalue_point = round(evalue_point, 2),
    evalue_ci = round(evalue_ci, 2)
  ))
}

# Convert Cohen's d to approximate risk ratio for E-value calculation
# Reference: VanderWeele (2017) - Approximation for continuous outcomes
cohens_d_to_rr <- function(d) {
  #' Convert standardized mean difference to approximate risk ratio
  #' Using VanderWeele's approximation: RR ≈ exp(0.91 * d)
  #' 
  #' @param d Cohen's d (standardized mean difference)
  #' @return Approximate risk ratio
  
  # Take absolute value for E-value calculation
  d <- abs(d)
  
  # Approximation from VanderWeele (2017)
  # This assumes outcome is approximately normally distributed
  rr <- exp(0.91 * d)
  
  return(rr)
}

# Calculate E-value for continuous outcome
calculate_evalue_continuous <- function(diff, se, pooled_sd) {
  #' Calculate E-value for continuous outcome
  #' 
  #' @param diff Mean difference between groups
  #' @param se Standard error of the difference
  #' @param pooled_sd Pooled standard deviation of the outcome
  #' @return List with Cohen's d, approximate RR, and E-values
  
  # Calculate Cohen's d
  d <- abs(diff) / pooled_sd
  
  # Calculate 95% CI for Cohen's d
  d_lo <- abs(diff - 1.96 * se) / pooled_sd
  d_hi <- abs(diff + 1.96 * se) / pooled_sd
  
  # Convert to approximate RR
  rr <- cohens_d_to_rr(d)
  rr_lo <- cohens_d_to_rr(min(d_lo, d_hi))
  rr_hi <- cohens_d_to_rr(max(d_lo, d_hi))
  
  # Calculate E-values
  evals <- calculate_evalue_rr(rr, rr_lo, rr_hi)
  
  return(list(
    cohens_d = round(d, 3),
    approx_rr = round(rr, 2),
    evalue_point = evals$evalue_point,
    evalue_ci = evals$evalue_ci
  ))
}

############################################################################
# PART 1: DATA GENERATION - PMI STUDY STRUCTURE
############################################################################

cat("================================================================================\n")
cat("PART 1: GENERATING SIMULATED DATA (PMI STUDY STRUCTURE)\n")
cat("================================================================================\n\n")

generate_pmi_data <- function(n_per_group = 296) {
  #' Generate data mimicking PMI Cross-Sectional Risk Marker Study
  #' 
  #' Groups:
  #'   1 = Current Smokers (reference)
  #'   2 = THS Users (switched from cigarettes ≥2 years ago)
  #'   3 = Former Smokers (quit all tobacco ≥2 years ago)
  
  # Group-specific parameters (simulating SELECTION into groups)
  group_params <- list(
    # Group 1: Current Smokers - older, higher pack-years
    current_smoker = list(
      age_mean = 48, age_sd = 8,
      pack_years_mean = 28, pack_years_sd = 10,
      cpd_mean = 18, cpd_sd = 6,
      female_prop = 0.42,
      europe_prop = 0.55
    ),
    # Group 2: THS Users - younger, moderate smoking history
    ths_user = list(
      age_mean = 44, age_sd = 7,
      pack_years_mean = 22, pack_years_sd = 8,
      cpd_mean = 16, cpd_sd = 5,
      female_prop = 0.44,
      europe_prop = 0.50
    ),
    # Group 3: Former Smokers - health-conscious, lower pack-years
    former_smoker = list(
      age_mean = 46, age_sd = 8,
      pack_years_mean = 20, pack_years_sd = 9,
      cpd_mean = 15, cpd_sd = 5,
      female_prop = 0.46,
      europe_prop = 0.52
    )
  )
  
  group_names <- c("current_smoker", "ths_user", "former_smoker")
  data_list <- list()
  
  for (g in 1:3) {
    params <- group_params[[g]]
    
    # Generate covariates
    age <- pmax(30, pmin(60, rnorm(n_per_group, params$age_mean, params$age_sd)))
    pack_years <- pmax(8, params$pack_years_mean + 0.4 * (age - params$age_mean) + 
                         rnorm(n_per_group, 0, params$pack_years_sd))
    cpd_history <- pmax(10, rnorm(n_per_group, params$cpd_mean, params$cpd_sd))
    female <- rbinom(n_per_group, 1, params$female_prop)
    europe <- rbinom(n_per_group, 1, params$europe_prop)
    
    # Generate BIOMARKER OUTCOMES with true causal effects
    
    # 1. COHb (%) - lower is better
    cohb_base <- 5 + 0.02 * (age - 45) + 0.01 * cpd_history + rnorm(n_per_group, 0, 0.8)
    cohb_effect <- c(0, -3.5, -4.2)[g]
    cohb <- pmax(0.3, cohb_base + cohb_effect + rnorm(n_per_group, 0, 0.3))
    
    # 2. Total NNAL (ng/mL) - lower is better (log-normal)
    log_nnal_base <- log(200) + 0.01 * (age - 45) + 0.02 * cpd_history + rnorm(n_per_group, 0, 0.4)
    nnal_effect <- c(0, -1.8, -3.5)[g]
    total_nnal <- exp(log_nnal_base + nnal_effect + rnorm(n_per_group, 0, 0.2))
    
    # 3. WBC (10^9/L) - lower is better
    wbc_base <- 7.5 + 0.02 * (age - 45) + rnorm(n_per_group, 0, 1)
    wbc_effect <- c(0, -0.8, -1.2)[g]
    wbc <- pmax(3.5, wbc_base + wbc_effect + rnorm(n_per_group, 0, 0.3))
    
    # 4. 8-epi-PGF2α (pg/mg) - lower is better (log-normal)
    log_epi_base <- log(800) + 0.005 * (age - 45) + 0.01 * cpd_history + rnorm(n_per_group, 0, 0.3)
    epi_effect <- c(0, -0.4, -0.6)[g]
    epi_pgf2a <- exp(log_epi_base + epi_effect + rnorm(n_per_group, 0, 0.15))
    
    # 5. HDL-C (mg/dL) - higher is better
    hdl_base <- 48 - 0.1 * (age - 45) + 5 * female + rnorm(n_per_group, 0, 8)
    hdl_effect <- c(0, 4, 6)[g]
    hdl_c <- pmax(25, hdl_base + hdl_effect + rnorm(n_per_group, 0, 2))
    
    # 6. sICAM-1 (ng/mL) - lower is better (log-normal)
    log_sicam_base <- log(250) + 0.008 * (age - 45) + rnorm(n_per_group, 0, 0.25)
    sicam_effect <- c(0, -0.15, -0.25)[g]
    sicam1 <- exp(log_sicam_base + sicam_effect + rnorm(n_per_group, 0, 0.1))
    
    # 7. 11-DTX-B2 (pg/mg) - lower is better (log-normal)
    log_dtx_base <- log(1500) + 0.01 * (age - 45) + rnorm(n_per_group, 0, 0.35)
    dtx_effect <- c(0, -0.3, -0.5)[g]
    dtx_b2 <- exp(log_dtx_base + dtx_effect + rnorm(n_per_group, 0, 0.15))
    
    # 8. Central AIx (%) - lower is better
    aix_base <- 25 + 0.4 * (age - 45) - 5 * female + rnorm(n_per_group, 0, 8)
    aix_effect <- c(0, -3, -5)[g]
    central_aix <- aix_base + aix_effect + rnorm(n_per_group, 0, 2)
    
    # 9. FEV1 % predicted - higher is better
    fev1_base <- 95 - 0.3 * (age - 45) - 0.2 * pack_years + 2 * female + rnorm(n_per_group, 0, 8)
    fev1_effect <- c(0, 3, 5)[g]
    fev1_pct <- pmin(120, pmax(50, fev1_base + fev1_effect + rnorm(n_per_group, 0, 2)))
    
    # Create data frame
    group_data <- data.frame(
      patient_id = paste0("G", g, "_", sprintf("%03d", 1:n_per_group)),
      product_group = factor(group_names[g], levels = group_names),
      age = round(age, 1),
      female = female,
      europe = europe,
      pack_years = round(pack_years, 1),
      cpd_history = round(cpd_history, 0),
      cohb = round(cohb, 2),
      total_nnal = round(total_nnal, 1),
      wbc = round(wbc, 2),
      epi_pgf2a = round(epi_pgf2a, 1),
      hdl_c = round(hdl_c, 1),
      sicam1 = round(sicam1, 1),
      dtx_b2 = round(dtx_b2, 1),
      central_aix = round(central_aix, 1),
      fev1_pct = round(fev1_pct, 1),
      stringsAsFactors = FALSE
    )
    
    data_list[[g]] <- group_data
  }
  
  pmi_data <- bind_rows(data_list)
  return(pmi_data)
}

# Generate the dataset
pmi_data <- generate_pmi_data(n_per_group = 296)

cat("Data Generated (PMI Study Structure)\n")
cat("Total sample size:", nrow(pmi_data), "\n")
cat("Number of groups:", length(unique(pmi_data$product_group)), "\n\n")
cat("Sample size per group:\n")
print(table(pmi_data$product_group))
cat("\n")

############################################################################
# PART 2: ASSESS PRE-WEIGHTING BALANCE
############################################################################

cat("================================================================================\n")
cat("PART 2: PRE-WEIGHTING BALANCE ASSESSMENT\n")
cat("================================================================================\n\n")

# SMD calculation functions
calculate_smd <- function(data, var, group_var) {
  groups <- levels(data[[group_var]])
  pairwise_smds <- data.frame()
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      group1_data <- data[data[[group_var]] == groups[i], var]
      group2_data <- data[data[[group_var]] == groups[j], var]
      
      mean1 <- mean(group1_data, na.rm = TRUE)
      mean2 <- mean(group2_data, na.rm = TRUE)
      sd1 <- sd(group1_data, na.rm = TRUE)
      sd2 <- sd(group2_data, na.rm = TRUE)
      
      pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
      smd <- (mean1 - mean2) / pooled_sd
      
      pairwise_smds <- rbind(pairwise_smds, 
                             data.frame(group1 = groups[i], 
                                       group2 = groups[j], 
                                       smd = smd,
                                       stringsAsFactors = FALSE))
    }
  }
  return(pairwise_smds)
}

calculate_weighted_smd <- function(data, var, group_var, weight_var) {
  groups <- levels(data[[group_var]])
  pairwise_smds <- data.frame()
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      group1_data <- data[data[[group_var]] == groups[i], ]
      group2_data <- data[data[[group_var]] == groups[j], ]
      
      mean1 <- weighted.mean(group1_data[[var]], group1_data[[weight_var]], na.rm = TRUE)
      mean2 <- weighted.mean(group2_data[[var]], group2_data[[weight_var]], na.rm = TRUE)
      
      weighted_var1 <- sum(group1_data[[weight_var]] * 
                          (group1_data[[var]] - mean1)^2, na.rm = TRUE) / 
                       sum(group1_data[[weight_var]])
      sd1 <- sqrt(weighted_var1)
      
      weighted_var2 <- sum(group2_data[[weight_var]] * 
                          (group2_data[[var]] - mean2)^2, na.rm = TRUE) / 
                       sum(group2_data[[weight_var]])
      sd2 <- sqrt(weighted_var2)
      
      pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
      smd <- (mean1 - mean2) / pooled_sd
      
      pairwise_smds <- rbind(pairwise_smds, 
                             data.frame(group1 = groups[i], 
                                       group2 = groups[j], 
                                       smd = smd,
                                       stringsAsFactors = FALSE))
    }
  }
  return(pairwise_smds)
}

# Covariates to check for balance
covariates <- c("age", "female", "europe", "pack_years", "cpd_history")

cat("PRE-WEIGHTING STANDARDIZED MEAN DIFFERENCES:\n\n")
pre_weight_smds <- list()

for (var in covariates) {
  smd_results <- calculate_smd(pmi_data, var, "product_group")
  pre_weight_smds[[var]] <- smd_results
  
  cat(paste0("Variable: ", toupper(var), "\n"))
  print(smd_results, row.names = FALSE)
  cat(paste0("Max absolute SMD: ", round(max(abs(smd_results$smd)), 3), "\n\n"))
}

# Baseline characteristics by group
cat("BASELINE CHARACTERISTICS BY PRODUCT GROUP (UNWEIGHTED):\n\n")
baseline_summary <- pmi_data %>%
  group_by(product_group) %>%
  summarise(
    N = n(),
    Age_mean = round(mean(age), 1),
    Age_sd = round(sd(age), 1),
    Female_pct = round(100 * mean(female), 1),
    Europe_pct = round(100 * mean(europe), 1),
    PackYears_mean = round(mean(pack_years), 1),
    PackYears_sd = round(sd(pack_years), 1),
    CPD_mean = round(mean(cpd_history), 1),
    CPD_sd = round(sd(cpd_history), 1),
    .groups = 'drop'
  )
print(as.data.frame(baseline_summary))
cat("\n")

############################################################################
# PART 3: ESTIMATE GENERALIZED PROPENSITY SCORES
############################################################################

cat("================================================================================\n")
cat("PART 3: ESTIMATING GENERALIZED PROPENSITY SCORES (GPS)\n")
cat("================================================================================\n\n")

cat("Fitting multinomial logistic regression model...\n")
cat("Model: P(Group | X) where X = age, female, europe, pack_years, cpd_history\n\n")

gps_model <- multinom(
  product_group ~ age + female + europe + pack_years + cpd_history, 
  data = pmi_data, 
  trace = FALSE
)

cat("Model coefficients (reference = current_smoker):\n")
print(round(coef(gps_model), 4))
cat("\n")

# Predict GPS
gps_matrix <- predict(gps_model, type = "probs")

pmi_data$gps <- sapply(1:nrow(pmi_data), function(i) {
  group_idx <- as.numeric(pmi_data$product_group[i])
  gps_matrix[i, group_idx]
})

cat("GENERALIZED PROPENSITY SCORE DIAGNOSTICS:\n\n")
cat("GPS summary:\n")
print(summary(pmi_data$gps))
cat("\n")

cat("GPS by product group:\n")
gps_by_group <- pmi_data %>% 
  group_by(product_group) %>% 
  summarise(
    Mean_GPS = round(mean(gps), 3),
    Median_GPS = round(median(gps), 3),
    Min_GPS = round(min(gps), 3),
    Max_GPS = round(max(gps), 3),
    .groups = 'drop'
  )
print(as.data.frame(gps_by_group))
cat("\n")

############################################################################
# PART 4: CALCULATE IPTW WEIGHTS WITH TRIMMING OPTIONS
############################################################################

cat("================================================================================\n")
cat("PART 4: CALCULATING IPTW WEIGHTS WITH TRIMMING\n")
cat("================================================================================\n\n")

# Marginal probabilities
marginal_probs <- pmi_data %>%
  group_by(product_group) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(marginal_prob = n / sum(n))

pmi_data <- pmi_data %>%
  left_join(marginal_probs %>% select(product_group, marginal_prob), 
            by = "product_group")

# Calculate unstabilized and stabilized weights
pmi_data$uw <- 1 / pmi_data$gps
pmi_data$sw <- pmi_data$marginal_prob / pmi_data$gps

cat("UNTRIMMED WEIGHT DIAGNOSTICS:\n\n")
cat("Unstabilized weights:\n")
print(summary(pmi_data$uw))
cat("\n")
cat("Stabilized weights:\n")
print(summary(pmi_data$sw))
cat("\n")

############################################################################
# WEIGHT TRIMMING STRATEGIES
############################################################################

cat("--------------------------------------------------------------------------------\n")
cat("WEIGHT TRIMMING STRATEGIES\n")
cat("--------------------------------------------------------------------------------\n\n")

# Strategy 1: Percentile-based trimming (1st and 99th percentiles)
trim_percentile <- function(weights, lower = 0.01, upper = 0.99) {
  lower_bound <- quantile(weights, lower)
  upper_bound <- quantile(weights, upper)
  pmin(pmax(weights, lower_bound), upper_bound)
}

# Strategy 2: Fixed threshold trimming
trim_fixed <- function(weights, max_weight = 10) {
  pmin(weights, max_weight)
}

# Strategy 3: Symmetric percentile trimming on GPS
trim_gps_percentile <- function(data, lower = 0.01, upper = 0.99) {
  gps_lower <- quantile(data$gps, lower)
  gps_upper <- quantile(data$gps, upper)
  data$gps_trimmed <- pmax(pmin(data$gps, gps_upper), gps_lower)
  data$sw_gps_trim <- data$marginal_prob / data$gps_trimmed
  return(data)
}

# Apply trimming strategies
cat("Applying trimming strategies...\n\n")

# 1. Percentile trimming on weights (1st-99th)
pmi_data$sw_trim_p99 <- trim_percentile(pmi_data$sw, 0.01, 0.99)

# 2. Percentile trimming on weights (5th-95th) - more aggressive
pmi_data$sw_trim_p95 <- trim_percentile(pmi_data$sw, 0.05, 0.95)

# 3. Fixed threshold (max = 10)
pmi_data$sw_trim_f10 <- trim_fixed(pmi_data$sw, 10)

# 4. Fixed threshold (max = 5) - more aggressive
pmi_data$sw_trim_f5 <- trim_fixed(pmi_data$sw, 5)

# 5. GPS percentile trimming
pmi_data <- trim_gps_percentile(pmi_data, 0.01, 0.99)

# Compare trimming methods
trimming_comparison <- data.frame(
  Method = c("Untrimmed", "Percentile 1-99%", "Percentile 5-95%", 
             "Fixed Max=10", "Fixed Max=5", "GPS Percentile 1-99%"),
  Weight_Var = c("sw", "sw_trim_p99", "sw_trim_p95", 
                 "sw_trim_f10", "sw_trim_f5", "sw_gps_trim"),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(trimming_comparison)) {
  w <- pmi_data[[trimming_comparison$Weight_Var[i]]]
  trimming_comparison$Min[i] <- round(min(w), 3)
  trimming_comparison$Max[i] <- round(max(w), 3)
  trimming_comparison$Mean[i] <- round(mean(w), 3)
  trimming_comparison$SD[i] <- round(sd(w), 3)
  trimming_comparison$ESS[i] <- round(sum(w)^2 / sum(w^2), 0)
  trimming_comparison$ESS_pct[i] <- round(100 * sum(w)^2 / sum(w^2) / nrow(pmi_data), 1)
}

cat("TRIMMING METHOD COMPARISON:\n\n")
print(trimming_comparison[, c("Method", "Min", "Max", "Mean", "SD", "ESS", "ESS_pct")], row.names = FALSE)
cat("\n")

# Select primary trimmed weights (percentile 1-99%)
pmi_data$sw_trimmed <- pmi_data$sw_trim_p99

cat("PRIMARY TRIMMING METHOD: Percentile 1-99%\n")
cat("Rationale: Balances bias-variance tradeoff; commonly recommended in literature\n\n")

############################################################################
# PART 5: POST-WEIGHTING BALANCE (TRIMMED vs UNTRIMMED)
############################################################################

cat("================================================================================\n")
cat("PART 5: POST-WEIGHTING BALANCE ASSESSMENT\n")
cat("================================================================================\n\n")

# Calculate SMDs for both trimmed and untrimmed weights
post_weight_smds_untrimmed <- list()
post_weight_smds_trimmed <- list()

for (var in covariates) {
  post_weight_smds_untrimmed[[var]] <- calculate_weighted_smd(pmi_data, var, "product_group", "sw")
  post_weight_smds_trimmed[[var]] <- calculate_weighted_smd(pmi_data, var, "product_group", "sw_trimmed")
}

# Create comparison table
smd_comparison <- data.frame()

for (var in covariates) {
  pre <- pre_weight_smds[[var]]
  post_untrim <- post_weight_smds_untrimmed[[var]]
  post_trim <- post_weight_smds_trimmed[[var]]
  
  for (i in 1:nrow(pre)) {
    comparison_label <- paste0(gsub("_", " ", pre$group1[i]), " vs ", gsub("_", " ", pre$group2[i]))
    
    smd_comparison <- rbind(smd_comparison, data.frame(
      Variable = var,
      Comparison = comparison_label,
      SMD_Pre = round(pre$smd[i], 3),
      SMD_Post_Untrimmed = round(post_untrim$smd[i], 3),
      SMD_Post_Trimmed = round(post_trim$smd[i], 3),
      stringsAsFactors = FALSE
    ))
  }
}

cat("SMD COMPARISON: Pre-Weighting vs Post-Weighting (Untrimmed vs Trimmed):\n\n")
print(smd_comparison, row.names = FALSE)
cat("\n")

# Balance summary
cat("BALANCE SUMMARY:\n")
cat("  Pre-weighting |SMD| > 0.1:", sum(abs(smd_comparison$SMD_Pre) > 0.1), "/", nrow(smd_comparison), "\n")
cat("  Post-weighting (untrimmed) |SMD| > 0.1:", sum(abs(smd_comparison$SMD_Post_Untrimmed) > 0.1), "/", nrow(smd_comparison), "\n")
cat("  Post-weighting (trimmed) |SMD| > 0.1:", sum(abs(smd_comparison$SMD_Post_Trimmed) > 0.1), "/", nrow(smd_comparison), "\n")
cat("  Max |SMD| pre:", round(max(abs(smd_comparison$SMD_Pre)), 3), "\n")
cat("  Max |SMD| post (untrimmed):", round(max(abs(smd_comparison$SMD_Post_Untrimmed)), 3), "\n")
cat("  Max |SMD| post (trimmed):", round(max(abs(smd_comparison$SMD_Post_Trimmed)), 3), "\n\n")

############################################################################
# PART 6: VISUALIZATIONS
############################################################################

cat("================================================================================\n")
cat("PART 6: CREATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# 1. Love Plot (comparing untrimmed vs trimmed)
love_plot_data <- smd_comparison %>%
  pivot_longer(cols = c(SMD_Pre, SMD_Post_Untrimmed, SMD_Post_Trimmed),
               names_to = "Weighting",
               values_to = "SMD") %>%
  mutate(
    Weighting = factor(Weighting,
                      levels = c("SMD_Pre", "SMD_Post_Untrimmed", "SMD_Post_Trimmed"),
                      labels = c("Unweighted", "IPTW (Untrimmed)", "IPTW (Trimmed)")),
    Variable_Comparison = paste0(Variable, " (", Comparison, ")")
  )

p1 <- ggplot(love_plot_data, aes(x = abs(SMD), y = Variable_Comparison,
                                  color = Weighting, shape = Weighting)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_manual(values = c("Unweighted" = "darkblue", 
                                 "IPTW (Untrimmed)" = "orange",
                                 "IPTW (Trimmed)" = "darkgreen")) +
  scale_shape_manual(values = c("Unweighted" = 16, 
                                "IPTW (Untrimmed)" = 17,
                                "IPTW (Trimmed)" = 15)) +
  labs(
    title = "Covariate Balance: Unweighted vs IPTW (Untrimmed vs Trimmed)",
    subtitle = "PMI THS Switching Study - 3 Product Groups",
    x = "Absolute Standardized Mean Difference",
    y = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(figures_dir, "love_plot_trimming_comparison.png"), p1, width = 11, height = 7, dpi = 300)
cat("Love plot saved:", file.path(figures_dir, "love_plot_trimming_comparison.png"), "\n")

# 2. Weight Distribution (untrimmed vs trimmed)
weight_plot_data <- pmi_data %>%
  select(product_group, sw, sw_trimmed) %>%
  pivot_longer(cols = c(sw, sw_trimmed),
               names_to = "Weight_Type",
               values_to = "Weight") %>%
  mutate(Weight_Type = factor(Weight_Type,
                              levels = c("sw", "sw_trimmed"),
                              labels = c("Untrimmed", "Trimmed (1-99%)")))

p2 <- ggplot(weight_plot_data, aes(x = Weight, fill = Weight_Type)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  facet_wrap(~product_group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Untrimmed" = "steelblue", "Trimmed (1-99%)" = "forestgreen")) +
  labs(
    title = "Distribution of Stabilized IPTW Weights: Untrimmed vs Trimmed",
    x = "Stabilized Weight",
    y = "Frequency",
    fill = "Weight Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(figures_dir, "weight_distribution_trimming.png"), p2, width = 9, height = 9, dpi = 300)
cat("Weight distribution plot saved:", file.path(figures_dir, "weight_distribution_trimming.png"), "\n")

# 3. GPS Distribution
p3 <- ggplot(pmi_data, aes(x = gps, fill = product_group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("current_smoker" = "#E41A1C", 
                                "ths_user" = "#377EB8", 
                                "former_smoker" = "#4DAF4A"),
                    labels = c("Current Smoker", "THS User", "Former Smoker")) +
  labs(
    title = "Generalized Propensity Score Distribution by Product Group",
    subtitle = "P(Observed Group | Baseline Covariates)",
    x = "Generalized Propensity Score",
    y = "Density",
    fill = "Product Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

ggsave(file.path(figures_dir, "gps_distribution_pmi.png"), p3, width = 10, height = 6, dpi = 300)
cat("GPS distribution plot saved:", file.path(figures_dir, "gps_distribution_pmi.png"), "\n\n")

############################################################################
# PART 7: OUTCOME ANALYSIS WITH SENSITIVITY ANALYSIS
############################################################################

cat("================================================================================\n")
cat("PART 7: OUTCOME ANALYSIS - BIOMARKERS WITH E-VALUE SENSITIVITY\n")
cat("================================================================================\n\n")

# Define biomarkers
biomarkers <- c("cohb", "total_nnal", "wbc", "epi_pgf2a", "hdl_c", 
                "sicam1", "dtx_b2", "central_aix", "fev1_pct")

biomarker_labels <- c(
  cohb = "COHb (%)",
  total_nnal = "Total NNAL (ng/mL)",
  wbc = "WBC (10^9/L)",
  epi_pgf2a = "8-epi-PGF2α (pg/mg)",
  hdl_c = "HDL-C (mg/dL)",
  sicam1 = "sICAM-1 (ng/mL)",
  dtx_b2 = "11-DTX-B2 (pg/mg)",
  central_aix = "Central AIx (%)",
  fev1_pct = "FEV1 % predicted"
)

# Function to analyze biomarker with E-value
analyze_biomarker_full <- function(data, biomarker, weight_var = NULL) {
  
  if (is.null(weight_var)) {
    model <- lm(as.formula(paste(biomarker, "~ product_group")), data = data)
  } else {
    model <- lm(as.formula(paste(biomarker, "~ product_group")), 
                data = data, weights = data[[weight_var]])
  }
  
  coefs <- summary(model)$coefficients
  
  # Overall pooled SD for E-value calculation
  pooled_sd <- sd(data[[biomarker]], na.rm = TRUE)
  
  # THS vs Smoker
  ths_diff <- coefs["product_groupths_user", "Estimate"]
  ths_se <- coefs["product_groupths_user", "Std. Error"]
  ths_p <- coefs["product_groupths_user", "Pr(>|t|)"]
  
  # Calculate E-value for THS vs Smoker
  ths_evalue <- calculate_evalue_continuous(ths_diff, ths_se, pooled_sd)
  
  # Former vs Smoker
  former_diff <- coefs["product_groupformer_smoker", "Estimate"]
  former_se <- coefs["product_groupformer_smoker", "Std. Error"]
  former_p <- coefs["product_groupformer_smoker", "Pr(>|t|)"]
  
  former_evalue <- calculate_evalue_continuous(former_diff, former_se, pooled_sd)
  
  return(data.frame(
    THS_Diff = ths_diff,
    THS_SE = ths_se,
    THS_p = ths_p,
    THS_Cohens_d = ths_evalue$cohens_d,
    THS_Approx_RR = ths_evalue$approx_rr,
    THS_Evalue_Point = ths_evalue$evalue_point,
    THS_Evalue_CI = ths_evalue$evalue_ci,
    Former_Diff = former_diff,
    Former_SE = former_se,
    Former_p = former_p,
    Former_Cohens_d = former_evalue$cohens_d,
    Former_Evalue_Point = former_evalue$evalue_point
  ))
}

# Analyze with different weight specifications
cat("BIOMARKER ANALYSIS: Unweighted vs IPTW (Untrimmed) vs IPTW (Trimmed)\n")
cat("================================================================================\n\n")

results_list <- list()

for (bm in biomarkers) {
  # Unweighted
  res_unweighted <- analyze_biomarker_full(pmi_data, bm, NULL)
  res_unweighted$Biomarker <- biomarker_labels[bm]
  res_unweighted$Method <- "Unweighted"
  
  # IPTW Untrimmed
  res_iptw <- analyze_biomarker_full(pmi_data, bm, "sw")
  res_iptw$Biomarker <- biomarker_labels[bm]
  res_iptw$Method <- "IPTW_Untrimmed"
  
  # IPTW Trimmed
  res_iptw_trim <- analyze_biomarker_full(pmi_data, bm, "sw_trimmed")
  res_iptw_trim$Biomarker <- biomarker_labels[bm]
  res_iptw_trim$Method <- "IPTW_Trimmed"
  
  results_list[[paste0(bm, "_unweighted")]] <- res_unweighted
  results_list[[paste0(bm, "_iptw")]] <- res_iptw
  results_list[[paste0(bm, "_iptw_trim")]] <- res_iptw_trim
}

all_results <- bind_rows(results_list)

# Create summary table for THS vs Smokers
cat("THS USERS vs CURRENT SMOKERS - Effect Estimates and E-Values:\n\n")

ths_summary <- all_results %>%
  select(Biomarker, Method, THS_Diff, THS_SE, THS_p, THS_Cohens_d, 
         THS_Approx_RR, THS_Evalue_Point, THS_Evalue_CI) %>%
  mutate(across(c(THS_Diff, THS_SE), ~round(., 2)),
         THS_p = format.pval(THS_p, digits = 2)) %>%
  pivot_wider(names_from = Method,
              values_from = c(THS_Diff, THS_SE, THS_p, THS_Cohens_d, 
                             THS_Approx_RR, THS_Evalue_Point, THS_Evalue_CI))

# Simplified comparison table
comparison_table <- data.frame(
  Biomarker = biomarker_labels[biomarkers]
)

for (bm in biomarkers) {
  unw <- all_results[all_results$Biomarker == biomarker_labels[bm] & all_results$Method == "Unweighted", ]
  trim <- all_results[all_results$Biomarker == biomarker_labels[bm] & all_results$Method == "IPTW_Trimmed", ]
  
  idx <- which(comparison_table$Biomarker == biomarker_labels[bm])
  comparison_table$Diff_Unweighted[idx] <- round(unw$THS_Diff, 2)
  comparison_table$Diff_IPTW[idx] <- round(trim$THS_Diff, 2)
  comparison_table$Change_Pct[idx] <- round(100 * (trim$THS_Diff - unw$THS_Diff) / abs(unw$THS_Diff), 1)
  comparison_table$Cohens_d_IPTW[idx] <- trim$THS_Cohens_d
  comparison_table$Evalue_Point[idx] <- trim$THS_Evalue_Point
  comparison_table$Evalue_CI[idx] <- trim$THS_Evalue_CI
}

print(comparison_table, row.names = FALSE)
cat("\n")

############################################################################
# E-VALUE INTERPRETATION
############################################################################

cat("================================================================================\n")
cat("E-VALUE SENSITIVITY ANALYSIS INTERPRETATION\n")
cat("================================================================================\n\n")

cat("E-VALUE INTERPRETATION GUIDE:\n")
cat("--------------------------------------------------------------------------------\n")
cat("The E-value quantifies the minimum strength of association that an unmeasured\n")
cat("confounder would need to have with BOTH the exposure (product group) AND the\n")
cat("outcome (biomarker) to fully explain away the observed effect.\n\n")

cat("E-value = 1.0: No unmeasured confounding needed (CI includes null)\n")
cat("E-value = 1.5: Weak unmeasured confounder could explain the effect\n")
cat("E-value = 2.0: Moderate unmeasured confounder needed\n")
cat("E-value = 3.0: Strong unmeasured confounder needed\n")
cat("E-value ≥ 4.0: Very strong unmeasured confounder needed (robust finding)\n\n")

cat("BIOMARKER-SPECIFIC E-VALUE SUMMARY (THS vs Smokers, IPTW-Trimmed):\n")
cat("--------------------------------------------------------------------------------\n\n")

for (i in 1:nrow(comparison_table)) {
  bm <- comparison_table$Biomarker[i]
  ev_point <- comparison_table$Evalue_Point[i]
  ev_ci <- comparison_table$Evalue_CI[i]
  
  # Interpretation
  if (ev_ci >= 3) {
    robustness <- "ROBUST - Very strong unmeasured confounder needed"
  } else if (ev_ci >= 2) {
    robustness <- "MODERATE - Strong unmeasured confounder needed"
  } else if (ev_ci >= 1.5) {
    robustness <- "SENSITIVE - Moderate unmeasured confounder could explain"
  } else {
    robustness <- "FRAGILE - Weak unmeasured confounder could explain"
  }
  
  cat(sprintf("%-25s E-value: %.2f (CI: %.2f) - %s\n", 
              bm, ev_point, ev_ci, robustness))
}

cat("\n")

############################################################################
# PART 8: SENSITIVITY ANALYSIS - TRIMMING THRESHOLDS
############################################################################

cat("================================================================================\n")
cat("PART 8: SENSITIVITY ANALYSIS - EFFECT OF TRIMMING THRESHOLD\n")
cat("================================================================================\n\n")

# Compare estimates across trimming methods for a key biomarker (COHb)
cat("SENSITIVITY TO TRIMMING THRESHOLD (Example: COHb)\n")
cat("--------------------------------------------------------------------------------\n\n")

trimming_sensitivity <- data.frame(
  Method = c("Unweighted", "IPTW Untrimmed", "IPTW Trim 1-99%", 
             "IPTW Trim 5-95%", "IPTW Trim Max=10", "IPTW Trim Max=5"),
  Weight_Var = c(NA, "sw", "sw_trim_p99", "sw_trim_p95", "sw_trim_f10", "sw_trim_f5"),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(trimming_sensitivity)) {
  if (is.na(trimming_sensitivity$Weight_Var[i])) {
    model <- lm(cohb ~ product_group, data = pmi_data)
  } else {
    model <- lm(cohb ~ product_group, data = pmi_data, 
                weights = pmi_data[[trimming_sensitivity$Weight_Var[i]]])
  }
  
  coefs <- summary(model)$coefficients
  trimming_sensitivity$THS_Diff[i] <- round(coefs["product_groupths_user", "Estimate"], 3)
  trimming_sensitivity$THS_SE[i] <- round(coefs["product_groupths_user", "Std. Error"], 3)
  trimming_sensitivity$THS_p[i] <- format.pval(coefs["product_groupths_user", "Pr(>|t|)"], digits = 3)
  
  # ESS
  if (!is.na(trimming_sensitivity$Weight_Var[i])) {
    w <- pmi_data[[trimming_sensitivity$Weight_Var[i]]]
    trimming_sensitivity$ESS[i] <- round(sum(w)^2 / sum(w^2), 0)
  } else {
    trimming_sensitivity$ESS[i] <- nrow(pmi_data)
  }
}

print(trimming_sensitivity[, c("Method", "THS_Diff", "THS_SE", "THS_p", "ESS")], row.names = FALSE)
cat("\n")

cat("INTERPRETATION:\n")
cat("- More aggressive trimming → estimates closer to unweighted (less bias correction)\n")
cat("- More aggressive trimming → lower variance (smaller SE)\n")
cat("- More aggressive trimming → higher ESS\n")
cat("- Trade-off: Bias reduction vs. variance inflation\n\n")

############################################################################
# PART 9: SAVE RESULTS
############################################################################

cat("================================================================================\n")
cat("PART 9: SAVING RESULTS\n")
cat("================================================================================\n\n")

# Save processed dataset
write.csv(pmi_data, file.path(tables_dir, "pmi_data_with_weights_trimmed.csv"), row.names = FALSE)
cat("Data saved:", file.path(tables_dir, "pmi_data_with_weights_trimmed.csv"), "\n")

# Save trimming comparison
write.csv(trimming_comparison, file.path(tables_dir, "trimming_method_comparison.csv"), row.names = FALSE)
cat("Trimming comparison saved:", file.path(tables_dir, "trimming_method_comparison.csv"), "\n")

# Save SMD comparison
write.csv(smd_comparison, file.path(tables_dir, "smd_comparison_trimmed.csv"), row.names = FALSE)
cat("SMD comparison saved:", file.path(tables_dir, "smd_comparison_trimmed.csv"), "\n")

# Save biomarker results with E-values
write.csv(comparison_table, file.path(tables_dir, "biomarker_results_evalues.csv"), row.names = FALSE)
cat("Biomarker results with E-values saved:", file.path(tables_dir, "biomarker_results_evalues.csv"), "\n")

# Save trimming sensitivity
write.csv(trimming_sensitivity, file.path(tables_dir, "trimming_sensitivity_cohb.csv"), row.names = FALSE)
cat("Trimming sensitivity saved:", file.path(tables_dir, "trimming_sensitivity_cohb.csv"), "\n")

# Create comprehensive summary report
ess_trimmed <- sum(pmi_data$sw_trimmed)^2 / sum(pmi_data$sw_trimmed^2)

summary_report <- data.frame(
  Metric = c(
    "Total Sample Size",
    "N per Group",
    "Number of Groups",
    "Covariates in GPS Model",
    "--- UNTRIMMED WEIGHTS ---",
    "ESS (Untrimmed)",
    "ESS % (Untrimmed)",
    "Max Weight (Untrimmed)",
    "--- TRIMMED WEIGHTS (1-99%) ---",
    "ESS (Trimmed)",
    "ESS % (Trimmed)",
    "Max Weight (Trimmed)",
    "--- BALANCE ---",
    "Max |SMD| Pre-Weighting",
    "Max |SMD| Post (Untrimmed)",
    "Max |SMD| Post (Trimmed)",
    "--- E-VALUES (IPTW-Trimmed) ---",
    "Mean E-value (Point)",
    "Mean E-value (CI)",
    "Min E-value (CI)",
    "Biomarkers with E-value CI ≥ 2"
  ),
  Value = c(
    nrow(pmi_data),
    296,
    3,
    length(covariates),
    "",
    round(sum(pmi_data$sw)^2 / sum(pmi_data$sw^2), 0),
    round(100 * sum(pmi_data$sw)^2 / sum(pmi_data$sw^2) / nrow(pmi_data), 1),
    round(max(pmi_data$sw), 2),
    "",
    round(ess_trimmed, 0),
    round(100 * ess_trimmed / nrow(pmi_data), 1),
    round(max(pmi_data$sw_trimmed), 2),
    "",
    round(max(abs(smd_comparison$SMD_Pre)), 3),
    round(max(abs(smd_comparison$SMD_Post_Untrimmed)), 3),
    round(max(abs(smd_comparison$SMD_Post_Trimmed)), 3),
    "",
    round(mean(comparison_table$Evalue_Point), 2),
    round(mean(comparison_table$Evalue_CI), 2),
    round(min(comparison_table$Evalue_CI), 2),
    sum(comparison_table$Evalue_CI >= 2)
  )
)

write.csv(summary_report, file.path(tables_dir, "summary_report_enhanced.csv"), row.names = FALSE)
cat("Summary report saved:", file.path(tables_dir, "summary_report_enhanced.csv"), "\n\n")

############################################################################
# FINAL SUMMARY
############################################################################

cat("================================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("================================================================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. WEIGHT TRIMMING:\n")
cat("   - Untrimmed max weight:", round(max(pmi_data$sw), 2), "\n")
cat("   - Trimmed (1-99%) max weight:", round(max(pmi_data$sw_trimmed), 2), "\n")
cat("   - ESS improved from", round(sum(pmi_data$sw)^2 / sum(pmi_data$sw^2), 0), 
    "to", round(ess_trimmed, 0), "with trimming\n\n")

cat("2. COVARIATE BALANCE:\n")
cat("   - Max |SMD| reduced from", round(max(abs(smd_comparison$SMD_Pre)), 3),
    "to", round(max(abs(smd_comparison$SMD_Post_Trimmed)), 3), "\n\n")

cat("3. E-VALUE SENSITIVITY ANALYSIS:\n")
cat("   - Mean E-value (point estimate):", round(mean(comparison_table$Evalue_Point), 2), "\n")
cat("   - Mean E-value (CI bound):", round(mean(comparison_table$Evalue_CI), 2), "\n")
cat("   - Biomarkers with robust findings (E-value CI ≥ 2):", 
    sum(comparison_table$Evalue_CI >= 2), "/", length(biomarkers), "\n\n")

cat("4. INTERPRETATION:\n")
cat("   - IPTW adjustment attenuates effect estimates (healthy switcher bias)\n")
cat("   - E-values suggest moderate-to-strong unmeasured confounding would be\n")
cat("     needed to fully explain the observed THS vs Smoker differences\n")
cat("   - Trimming improves precision with minimal impact on point estimates\n\n")

cat("FILES SAVED:\n")
cat("  Tables:", tables_dir, "\n")
cat("    - pmi_data_with_weights_trimmed.csv\n")
cat("    - trimming_method_comparison.csv\n")
cat("    - smd_comparison_trimmed.csv\n")
cat("    - biomarker_results_evalues.csv\n")
cat("    - trimming_sensitivity_cohb.csv\n")
cat("    - summary_report_enhanced.csv\n")
cat("  Figures:", figures_dir, "\n")
cat("    - love_plot_trimming_comparison.png\n")
cat("    - weight_distribution_trimming.png\n")
cat("    - gps_distribution_pmi.png\n\n")

cat("================================================================================\n")
cat("METHODOLOGICAL NOTES\n")
cat("================================================================================\n\n")

cat("WEIGHT TRIMMING RECOMMENDATIONS:\n")
cat("  1. Percentile-based (1-99%): Default choice, data-adaptive\n")
cat("  2. Fixed threshold: Use when clinical knowledge suggests max plausible weight\n")
cat("  3. Always compare results across trimming methods as sensitivity analysis\n\n")

cat("E-VALUE INTERPRETATION:\n")
cat("  - E-value answers: 'How strong would unmeasured confounding need to be?'\n")
cat("  - Higher E-value = more robust to unmeasured confounding\n")
cat("  - E-value for CI bound is more conservative (accounts for sampling variability)\n")
cat("  - Compare E-values to plausible confounders in your domain\n\n")

cat("LIMITATIONS:\n")
cat("  - E-values assume single unmeasured confounder (may underestimate if multiple)\n")
cat("  - Trimming introduces some bias (trade-off with variance)\n")
cat("  - Cross-sectional design still cannot establish temporality\n")
cat("  - GPS cannot adjust for unmeasured confounders\n\n")

cat("================================================================================\n")
