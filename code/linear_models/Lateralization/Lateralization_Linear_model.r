# Description: This script loads essential R packages for conducting linear and generalized linear modeling,
library(mgcv) # Generalized Additive Models (GAMs) - for smooth age effects
library(lme4) # Linear Mixed-Effects Models - alternative approach
library(lmerTest) # P-values for lmer models
library(ggplot2) # High-quality graphics
library(dplyr) # Data manipulation
library(tidyr) # Data reshaping
library(emmeans) # Estimated marginal means for post-hoc tests
library(gratia) # GAM visualization and diagnostics
library(patchwork) # Combining multiple plots
library(car) # Variance inflation factors, Levene's test
library(moments) # Skewness and kurtosis
library(effectsize) # Effect size calculations (Cohen's d)
library(performance) # Model diagnostics and ICC
library(gt) # Publication-ready tables
library(knitr) # Alternative table formatting
library(scales) # For formatting plot axes
library(clubSandwich)
library(lmtest)
library(pwr)
library(geepack)
library(ICC)
library(ggeffects) # Model predictions


#Set random seed for reproducibility
set.seed(333)

# Create output directory
setwd("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/q1k_neurosubs/code/linear_models")
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)
dir.create("output/diagnostics", showWarnings = FALSE)

# Load Data
data_path <- "C:/Users/gabot/OneDrive - McGill University/Desktop/github_repos/q1k_neurosubs/outputs/eeg/eeg_data.csv"
raw_data <- read.csv(data_path, stringsAsFactors = FALSE)

print(head(raw_data))
# Inspect data structure

required_cols <- c("subject", "fam_id", "affected_group", "eeg_age", "sex", "site",
                   "l_index_frontal")
missing_cols <- setdiff(required_cols, colnames(raw_data))
if (length(missing_cols) > 0) {
  stop(paste("The following required columns are missing from the data:", paste(missing_cols, collapse = ", ")))
}


# Working dataset
data <- raw_data

# Subset only siblings
data <- data %>% 
  filter(!group %in% c("mother", "father"))

# Affected group factor levels
data$affected_group <- factor(
  data$affected_group,
  levels = c("non-affected", "affected", "asd"),
  labels = c("Non-Affected", "Affected (non-ASD)", "ASD")
)
print(table(data$affected_group))

# Sex factor levels
data$sex <- factor(data$sex, levels = c("female", "male"))
print(table(data$sex))

# Site factor levels
data$site <- factor(data$site, levels = c("hsj", "mni"))
print(table(data$site))

# Family ID as factor
data$fam_id <- factor(data$fam_id)
print(paste("Number of unique families: ", length(levels(data$fam_id))))  

#Ensure age is numeric 
data$eeg_age <- as.numeric(data$eeg_age)
print(summary(data$eeg_age))


# Ensure outcome variables are numeric
 outcome_vars <- c("l_index", "l_index_frontal","pwr_frontal_rel")
for (var in outcome_vars) {
  data[[var]] <- as.numeric(data[[var]])
}

 
### SECTION 1:  Missing Data Analysis
 missing_summary <- data.frame(
   Variable = colnames(data),
   N_Missing = sapply(data, function(x) sum(is.na(x))),
   Percent_Missing = sapply(data, function(x) round(100 * sum(is.na(x)) / length(x), 2))
 )
 rownames(missing_summary) <- NULL
 
 print("Missing data summary:")
 print(missing_summary)
 
### SECTION 2 Outlier detection
 
 detect_outliers <- function(x, threshold = 3) {
   if (all(is.na(x))) return(rep(FALSE, length(x)))
   mean_x <- mean(x, na.rm = TRUE)
   sd_x <- sd(x, na.rm = TRUE)
   abs(x - mean_x) > threshold * sd_x
 }
 
 outlier_summary <- data.frame(Variable = outcome_vars, N_Outliers = NA, Percent_Outliers = NA)
 
 for (i in 1:length(outcome_vars)) {
   var <- outcome_vars[i]
   outliers <- detect_outliers(data[[var]])
   data[[paste0(var, "_outlier")]] <- outliers
   
   n_outliers <- sum(outliers, na.rm = TRUE)
   pct_outliers <- round(100 * n_outliers / sum(!is.na(data[[var]])), 2)
   
   outlier_summary$N_Outliers[i] <- n_outliers
   outlier_summary$Percent_Outliers[i] <- pct_outliers
  
 }
 
 print("\nOutlier summary:")
 print(outlier_summary)
 
 # Remove outliers for log_psd
 data <- data[!data$l_index_frontal_outlier , ]
 
 # SECTION 3:  Descriptive statistics
 
 print("OVERALL SAMPLE CHARACTERISTICS:")
 print(paste("Total N:", nrow(data)))
 print(paste("Age: M =", round(mean(data$eeg_age, na.rm = TRUE), 1),
             "SD =", round(sd(data$eeg_age, na.rm = TRUE), 1),
             "Range:", round(min(data$eeg_age, na.rm = TRUE), 1), "-",
             round(max(data$eeg_age, na.rm = TRUE), 1), "years"))
 print(paste("Sex:", sum(data$sex == "male", na.rm = TRUE), "male,",
             sum(data$sex == "female", na.rm = TRUE), "female"))
 print(paste("Site:", sum(data$site == "hsj", na.rm = TRUE), "hsj,",
             sum(data$site == "mni", na.rm = TRUE), "mni"))
 
 get_descriptives <- function(x) {
   c(
     N = sum(!is.na(x)),
     Mean = mean(x, na.rm = TRUE),
     SD = sd(x, na.rm = TRUE),
     Min = min(x, na.rm = TRUE),
     Max = max(x, na.rm = TRUE),
     Median = median(x, na.rm = TRUE),
     Q25 = quantile(x, 0.25, na.rm = TRUE),
     Q75 = quantile(x, 0.75, na.rm = TRUE)
   )
 }
 
 
 for (group in levels(data$affected_group)) {
   print(paste("\n========== GROUP:", group, "=========="))
   
   group_data <- data[data$affected_group == group, ]
   
   print(paste("N =", nrow(group_data)))
   print(paste("Age: M =", round(mean(group_data$eeg_age, na.rm = TRUE), 1),
               "SD =", round(sd(group_data$eeg_age, na.rm = TRUE), 1)))
   print(paste("Sex:", sum(group_data$sex == "male"), "male,",
               sum(group_data$sex == "female"), "female"))
   print(paste("Age groups:",
               sum(group_data$developmental_group == "Child (5-12y)"), "child,",
               sum(group_data$developmental_group == "Adolescent (13-17y)"), "adolescent,",
               sum(group_data$developmental_group == "Adult (18+y)"), "adult"))
   
   print("\nOutcome variables:")
   desc_table <- t(sapply(outcome_vars, function(v) get_descriptives(group_data[[v]])))
   print(round(desc_table, 2))
 }
 
 # Comprehensive table
 desc_full <- data %>%
   group_by(affected_group) %>%
   summarise(
     N = n(),
     Age_Mean = mean(eeg_age, na.rm = TRUE),
     Age_SD = sd(eeg_age, na.rm = TRUE),
     Pct_Male = 100 * mean(sex == "male"),
  #   Alpha_Mean = mean(alpha_Power, na.rm = TRUE),
   #  Alpha_SD = sd(alpha_Power, na.rm = TRUE),
    # Theta_Mean = mean(theta_Power, na.rm = TRUE),
     #Theta_SD = sd(theta_Power, na.rm = TRUE),
    # Disengage_Mean = mean(disengagement, na.rm = TRUE),
    # Disengage_SD = sd(disengagement, na.rm = TRUE),
    # Facilit_Mean = mean(facilitation, na.rm = TRUE),
   #  Facilit_SD = sd(facilitation, na.rm = TRUE),
     .groups = "drop"
   )
write.csv(desc_full, "output/tables/descriptive_statistics_by_group.csv", row.names = FALSE)

# Check family patterns
family_structure <- data %>%
  group_by(fam_id) %>%
  summarise(
    n_members = n(),
    n_asd = sum(affected_group == "ASD"),
    n_affected = sum(affected_group == "Affected (non-ASD)"),
    n_unaffected = sum(affected_group == "Non-Affected"),
    age_range = max(eeg_age) - min(eeg_age),
    has_multiple = n() > 1
  )

# Summary statistics
cat("Family structure summary:\n")
cat("Total families:", nrow(family_structure), "\n")
cat("Families with 1 member:", sum(family_structure$n_members == 1), "\n")
cat("Families with 2+ members:", sum(family_structure$n_members >= 2), "\n")
cat("Median family size:", median(family_structure$n_members), "\n")


# SECTION 3.1: AGE DIFFERENCES BETWEEN GROUPS
print("TESTING AGE DIFFERENCES BETWEEN GROUPS")

age_anova <- aov(eeg_age ~ affected_group, data = data)
print(summary(age_anova))
# Save ANOVA test summary
capture.output(summary(age_anova), file = "output/tables/age_anova_summary.txt")

# Post-hoc tests
age_posthoc <- TukeyHSD(age_anova)
print("\nPairwise age comparisons (Tukey HSD):")
print(age_posthoc)

# Plot Age Differences
ggplot(data, aes(x=eeg_age, fill=affected_group)) +
  geom_density(alpha=0.5) +
  labs(title="Age Distributions Reveal Minimal Overlap")+ 
  theme_minimal()

# Save plot
ggsave("./output/figures/age_distribution_by_group.png",  bg = "white",
       width=8, height=6, dpi=300)

# Save post-hoc test summary
capture.output(age_posthoc, file = "output/tables/age_posthoc_summary.txt")

# SECTION 4: ASSUMPTION CHECKING - DISTRIBUTIONS

outcome <- "l_index_frontal"

data_complete <- data[!is.na(data[[outcome]]), ]
print("Overall distribution:")
print(paste(" Skewness:", round(skewness(data_complete[[outcome]]), 3)))
print(paste(" Kurtosis:", round(kurtosis(data_complete[[outcome]]), 3)))
shapiro_test <- shapiro.test(data_complete[[outcome]])
# Shapiro-Wilk test
print(paste(" Shapiro-Wilk test: W =", round(shapiro_test$statistic, 4),
            ", p =", format.pval(shapiro_test$p.value, digits = 3)))

# Kolmogorov-Smirnov test against normal distribution
ks_test <- ks.test(data_complete[[outcome]], "pnorm",
                   mean = mean(data_complete[[outcome]]),
                   sd = sd(data_complete[[outcome]]))
print(paste(" Kolmogorov-Smirnov test: D =", round(ks_test$statistic, 4),
            ", p =", format.pval(ks_test$p.value, digits = 3)))


# Normality by group 
print("\nNormality by group:")
for (grp in levels(data_complete$affected_group)) {
  grp_data <- data_complete[data_complete$affected_group == grp, ][[outcome]]
  if (length(grp_data) >= 3 && length(grp_data) < 5000) {
    shap_grp <- shapiro.test(grp_data)
    print(paste(" ", grp, ": W =", round(shap_grp$statistic, 4),
                ", p =", format.pval(shap_grp$p.value, digits = 3)))
  }}

# SECTION 5 : Diagnostic plots

p1 <- ggplot(data_complete, aes(x = .data[[outcome]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(title = paste(outcome, "- Overall Distribution"),
       x = outcome, y = "Density") +
  theme_minimal(base_size = 10)

p2 <- ggplot(data_complete, aes(sample = .data[[outcome]])) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot (Overall)", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 10)

p3 <- ggplot(data_complete, aes(sample = .data[[outcome]], color = affected_group)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~affected_group) +
  labs(title = "Q-Q Plots by Group", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

p4 <- ggplot(data_complete, aes(x = affected_group, y = .data[[outcome]], fill = affected_group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  labs(title = "Distribution by Group", x = "Group", y = outcome) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
# Combined plot
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_annotation(title = paste("Distribution Diagnostics:", outcome),
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))
combined_plot
# Save plot

ggsave(filename = paste0("output/diagnostics/distribution_", outcome, ".png"),
       bg = "white",
       plot = combined_plot, width = 12, height = 10, dpi = 300)

print(paste("Distribution diagnostic plot saved for", outcome))


# SECTION 6: Homogeneity of Variance

# Levenes test
levene_result <- car::leveneTest(
  as.formula(paste(outcome, "~ affected_group")),
  data = data_complete,
  center = median # More robust to outliers than mean
)

print(levene_result)

# Variance by group
var_by_group <- data_complete %>%
  group_by(affected_group) %>%
  summarise(
    Variance = var(.data[[outcome]], na.rm = TRUE),
    SD = sd(.data[[outcome]], na.rm = TRUE),
    .groups = "drop"
  )

print("\nVariance by group:")
print(var_by_group)

# Variance Ratio
var_ratio <- max(var_by_group$Variance) / min(var_by_group$Variance)
print(paste("\nVariance ratio (max/min):", round(var_ratio, 2)))
if (var_ratio > 4) {
  print("Warning: Variance ratio exceeds 4, indicating potential heteroscedasticity.")
} else {
  print("Variance ratio is within acceptable limits.")
}

# Save Levene's test result
results_file <- paste0("output/tables/levene_test_", outcome, ".txt")

p_var <- ggplot(data_complete, aes(x = affected_group, y = .data[[outcome]])) +
  geom_violin(aes(fill = affected_group), alpha = 0.3, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(aes(color = affected_group), width = 0.2, alpha = 0.4, size = 1) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red", shape = 18) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "red") +
  theme_minimal()+
  labs(title = paste("Variance Structure:", outcome),
       subtitle = paste("Levene's test p =", round(levene_result$`Pr(>F)`[1], 4)))
p_var
ggsave(filename = paste0("output/diagnostics/variance_", outcome, ".png"), bg = "white", 
       plot = p_var, width = 8, height = 6, dpi = 300)


# SECTION 7: AGE OUTCOME RELATIONSHIP


# New plot
p_raw <- ggplot(data, aes(x = eeg_age, y = .data[[outcome]], color = affected_group)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), 
              se = TRUE, linewidth = 1.2) +
  scale_color_manual(
    values = c("Non-Affected" = "#2ecc71",
               "Affected (non-ASD)" = "#3498db",
               "ASD" = "#e74c3c")
  ) +
  labs(
    subtitle = "LOESS curves show observed trajectories",
    x = "Age (years)",
    y = "Lateralization",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

p_raw

# Plot raw with restricted age
data_age_restricted <- data %>% 
  filter(eeg_age >= common_min & eeg_age <= 18)

p_raw_restricted <- ggplot(data_age_restricted, aes(x = eeg_age, y = data_age_restricted[[outcome]], color = affected_group)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), 
              se = TRUE, linewidth = 1.2) +
  labs(
    x = "Age (years)",
    y = "Lateralization",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

p_raw_restricted


# Check association with age for categorical variables 
print("\nAge distribution by sex:")
age_by_sex <- data %>%
  group_by(sex) %>%
  summarise(Mean_Age = mean(eeg_age), SD_Age = sd(eeg_age), .groups = "drop")
print(age_by_sex)

t_test_sex_age <- t.test(eeg_age ~ sex, data = data)
print(paste("t-test for age difference by sex: p =", round(t_test_sex_age$p.value, 4)))

print("\nAge distribution by site:")
age_by_site <- data %>%
  group_by(site) %>%
  summarise(Mean_Age = mean(eeg_age), SD_Age = sd(eeg_age), .groups = "drop")
print(age_by_site)

t_test_site_age <- t.test(eeg_age ~ site, data = data)
print(paste("t-test for age difference by site: p =", round(t_test_site_age$p.value, 4)))

# Check out sex differences in each group 
affected_group_by_sex <- table(data$affected_group, data$sex)
affected_group_by_sex

# SECTION 8: Linear ModelS (AGE)

# 8.0 NULL MOdel (ICC)

formula_m0 <- as.formula(paste(outcome, "~ 1 + (1|fam_id)"))
m0 <- lmer(formula_m0, data = data, REML = TRUE)

icc_m0 <- performance::icc(m0)
cat("ICC (Family):", round(icc_m0$ICC_adjusted, 3), "\n")
cat("Interpretation:", round(100 * icc_m0$ICC_adjusted, 1), 
    "% of variance is between families\n\n")

# Interpretation
# - ICC < 0.05: Weak family effects, random intercept may not be needed
# - ICC 0.05-0.15: Moderate family effects
# - ICC > 0.15: Strong family effects -


# Check for VIF

# Rememebrs that High VIF (>5-10) is problematic 
prelim_model <- lm(l_index_frontal ~ eeg_age + affected_group + sex + site, data = data_complete)
vif_values <- car::vif(prelim_model)
print(vif_values)
# Save VIF results
write.csv(vif_values, "output/tables/vif_results.csv", row.names = FALSE)

#8.1  Linear model
mod_linear <- lm(as.formula(paste(outcome, "~ eeg_age + affected_group")),
                 data = data_complete)
aic_linear <- AIC(mod_linear)
print(paste("Linear model AIC:", round(aic_linear, 2)))
# Summary
summary_linear <- summary(mod_linear)
print(summary_linear)


# 8.2 Quadratic model
mod_quad <- lm(as.formula(paste(outcome, "~ eeg_age + I(eeg_age^2) + affected_group")),
               data = data_complete)
aic_quad <- AIC(mod_quad)
print(paste("Quadratic model AIC:", round(aic_quad, 2)))
# Summary
summary_quad <- summary(mod_quad)
print(summary_quad)

# 8.3 GAM model
mod_gam <- gam(as.formula(paste(outcome, "~ s(eeg_age, k = 5) + affected_group")),
               data = data_complete)
aic_gam <- AIC(mod_gam)
print(paste("GAM model AIC:", round(aic_gam, 2)))
# Summary
summary_gam <- summary(mod_gam)
print(summary_gam)

# Comparisons
aic_comparison <- data.frame(
  Model = c("Linear", "Quadratic", "GAM"),
  AIC = c(aic_linear, aic_quad, aic_gam),
  Delta_AIC = c(aic_linear, aic_quad, aic_gam) - min(c(aic_linear, aic_quad, aic_gam))
)
print("\nAIC Comparison:")
print(aic_comparison)

# Save comparison table

best_model <- aic_comparison$Model[which.min(aic_comparison$AIC)]
delta_aic <- min(aic_comparison$Delta_AIC[aic_comparison$Model != best_model])

print("**BEST MODEL SEEMS TO BE GAM***")


# SECTION 10: ****MAIN ANALYSIS****

# Save summaries
model_results <- list()
model_summaries <- list()

# Use only complete cases
data_model <- data[!is.na(data[[outcome]]), ]
n_complete <- nrow(data_model)
print(paste("Complete cases for", outcome, ":", n_complete))

# Save outcome name
model_results[[outcome]] <- list()

# MODEL 1: NULL MODEL (random effects only)

## USE BELOW IF GAM IS BETTER, IF NOT USE LMER


formula_m1 <- as.formula(paste(outcome, "~ 1 + s(fam_id, bs='re')"))

m1 <- gam(formula_m1, data = data_model, method = "REML")

model_results[[outcome]]$m1 <- m1

# Print model summary
print(summary(m1))

#extract key statistics
aic_m1 <- AIC(m1)
bic_m1 <- BIC(m1)
dev_exp_m1 <- summary(m1)$dev.expl * 100

print(paste("AIC:", round(aic_m1, 2)))
print(paste("BIC:", round(bic_m1, 2)))
print(paste("Deviance Explained:", round(dev_exp_m1, 2), "%"))

# ICC = between group variance/total variance
re_var <- gam.vcomp(m1)

print(re_var)

# Row 1 = family random effect variance
# Row 2 = residual variance 

var_family <- re_var[1, 1]^2  # Square the std.dev to get variance
var_residual <- re_var[2, 1]^2

total_var <- var_family + var_residual
icc_family <- var_family / total_var

cat("Family variance:", round(var_family, 6), "\n")
cat("Residual variance:", round(var_residual, 6), "\n")
cat("ICC (proportion due to family):", round(icc_family, 4), "\n")
cat("=> ", round(icc_family * 100, 2), "% of variance is due to family clustering\n\n")

# BASED ON THIS REMOVE THE FAM_ID FROM MODEL

# MODEL 2 MAIN EFFECTS (No interaction)


formula_m2 <- as.formula(paste(outcome,
                               "~ affected_group + s(eeg_age, k=5) + sex + site"))

m2 <- gam(formula_m2, data = data_model, method = "REML")

model_results[[outcome]]$m2 <- m2
# Print model summary
print(summary(m2))

aic_m2 <- AIC(m2)
bic_m2 <- BIC(m2)
dev_exp_m2 <- summary(m2)$dev.expl * 100

print(paste("AIC:", round(aic_m2, 2)))
print(paste("BIC:", round(bic_m2, 2)))
print(paste("Deviance explained:", round(dev_exp_m2, 2), "%"))
print(paste("Improvement over null:", round(dev_exp_m2 - dev_exp_m1, 2), "percentage points"))

# MODEL 3 Age x Group interaction (MAIN MODEL)

formula_m3 <- as.formula(paste(outcome,
                               "~ affected_group + s(eeg_age, by=affected_group, k=5) + sex + site"))

m3 <- gam(formula_m3, data = data_model, method = "REML")

model_results[[outcome]]$m3 <- m3
print(summary(m3))

aic_m3 <- AIC(m3)
bic_m3 <- BIC(m3)
dev_exp_m3 <- summary(m3)$dev.expl * 100

print(paste("AIC:", round(aic_m3, 2)))
print(paste("BIC:", round(bic_m3, 2)))
print(paste("Deviance explained:", round(dev_exp_m3, 2), "%"))


# Model 4 SEX INTERACTION

formula_m4 <- as.formula(paste(outcome,
                               "~ affected_group * sex + s(eeg_age, by=interaction(affected_group,sex),
                               k=4) +site"))

m4 <- gam(formula_m4, data = data_model, method = "REML")

model_results[[outcome]]$m4 <- m4
print(summary(m4))

aic_m4 <- AIC(m4)
bic_m4 <- BIC(m4)
dev_exp_m4 <- summary(m4)$dev.expl * 100

print(paste("AIC:", round(aic_m4, 2)))
print(paste("BIC:", round(bic_m4, 2)))
print(paste("Deviance explained:", round(dev_exp_m4, 2), "%"))


# MODEL COMPARISON TABLE
comparison_table <- data.frame(
  Model = c("M1: Null", "M2: Main Effects", "M3: Interaction x Group (GAM)", "M4:Sex INteraction  (lmer)"),
  AIC = c(aic_m1, aic_m2, aic_m3, aic_m4),
  BIC = c(bic_m1, bic_m2, bic_m3, bic_m4),
  Dev_Explained_or_R2 = c(
    paste0(round(dev_exp_m1, 1), "%"),
    paste0(round(dev_exp_m2, 1), "%"),
    paste0(round(dev_exp_m3, 1), "%"),
    paste0(round(dev_exp_m4, 1), "%")
    
  )
)

comparison_table$Delta_AIC <- comparison_table$AIC - min(comparison_table$AIC)
comparison_table$Delta_BIC <- comparison_table$BIC - min(comparison_table$BIC)

print(comparison_table)

# Save table comparison 
model_summaries[[outcome]]$comparison <- comparison_table

# Formal LRT
cat("\n--- Likelihood Ratio Tests ---\n")
cat("M2 vs M3 (does age×group interaction help?):\n")
print(anova(m2, m3, test = "Chisq"))

cat("\nM3 vs M4 (does sex×group interaction help?):\n")
print(anova(m3, m4, test = "Chisq"))


# Select the best model 
model_results[[outcome]]$best <- best_model

# SECTION 12: BEST MODEL DIAGNOSTICS

# Extract best model
best_model <- m4

data_model <- data[!is.na(data[[outcome]]), ]

# OVerall performance 

# 12.1 Performance metrics

perf <- model_performance(best_model)
print("Overall Model Performance:")
print(perf)

# Save to table
perf_table <- as.data.frame(perf)
write.csv(perf_table, 
          paste0("output/tables/model_performance_", outcome, ".csv"), 
          row.names = FALSE)


# 12.2. R-SQUARED DECOMPOSITION
# Marginal R² = variance explained by FIXED effects only
# Conditional R² = variance explained by FIXED + RANDOM effects

r2_vals <- performance::r2(best_model)
r2_vals

resid_df <- data.frame(
  fitted = fitted(best_model),
  residuals = residuals(best_model),
  eeg_age = data$eeg_age,
  affected_group = data$affected_group,
  sex = data$sex,
  fam_id = data$fam_id
)


# A. Normality of residuals
shapiro_resid <- shapiro.test(resid_df$residuals)
cat("\nShapiro-Wilk test on residuals: W =", round(shapiro_resid$statistic, 4),
    ", p =", round(shapiro_resid$p.value, 4), "\n")

# B. Homoscedasticity
bp_test <- lmtest::bptest(lm(residuals ~ fitted, data = resid_df))
cat("Breusch-Pagan test: χ² =", round(bp_test$statistic, 2),
    ", p =", round(bp_test$p.value, 4), "\n")

# C. Influential observations
check_outliers(best_model)

# D. Multicollinearity (re-check with final model)
check_collinearity(best_model)

# E. Overall diagnostic plot
check_model(best_model, residual_type = "normal")
ggsave("output/diagnostics/model_diagnostics_complete.png",
       width = 12, height = 10, dpi = 300)


# G. Cook's distance for influential cases
cooksd <- cooks.distance(lm(residuals ~ fitted, data = resid_df))
influential <- which(cooksd > 4/nrow(resid_df))
cat("\nInfluential cases (Cook's D > 4/n):", length(influential), "\n")
if (length(influential) > 0) {
  cat("IDs:", head(influential, 10), "\n")
}

# Print cohen's D for Groups
group_means <- data_complete %>%
  group_by(affected_group) %>%
  summarise(M = mean(l_index_frontal), SD = sd(l_index_frontal))

cohen_d_asd_control <- (group_means$M[3] - group_means$M[1]) / 
  sqrt(mean(c(group_means$SD[1]^2, group_means$SD[3]^2)))
print(paste("Cohen's d:", round(cohen_d_asd_control, 3)))


# REPORT MAIN EFFECTS

emm_group <- emmeans(m4, ~ affected_group)

cat("\n=== ESTIMATED MARGINAL MEANS BY GROUP ===\n")
print(emm_group)

# Pairwise comparisons
pairs_group <- pairs(emm_group, adjust = "tukey")
cat("\n=== PAIRWISE COMPARISONS (Tukey-adjusted) ===\n")
print(pairs_group)

# Effect sizes
eff_sizes <- eff_size(emm_group, sigma = sigma(best_model), edf = df.residual(best_model))
cat("\n=== EFFECT SIZES (Cohen's d) ===\n")
print(eff_sizes)

# INTERPRETATION
# - Cohen's d quantifies magnitude: |d| < 0.2 = trivial, 0.2-0.5 = small, 
#   0.5-0.8 = medium, >0.8 = large



## SENSITIVITY ANALYSI 


# Identify overlapping age range
age_overlap <- data %>%
  group_by(affected_group) %>%
  summarise(min_age = min(eeg_age), max_age = max(eeg_age))
age_overlap

common_min <- max(age_overlap$min_age)
common_max <- min(age_overlap$max_age)

cat("Common age range:", common_min, "to", common_max, "\n")

# Subset data
data_age_restricted <- data %>% 
  filter(eeg_age >= common_min & eeg_age <= 18)

cat("Sample sizes after restriction:\n")
print(table(data_age_restricted$affected_group))

# Re-run Model on this subset

best_age_restricted <- gam(formula_m3,
  data = data_age_restricted, method = "REML"
)

  # Compare to full sample
summary(best_model)  # Original

summary(best_age_restricted)  # Age-restricted

# COMPARE
emm_original <- emmeans(best_model, ~ affected_group)
emm_restricted <- emmeans(best_age_restricted, ~ affected_group)

cat("\nEffect Sizes - Original vs Age-Restricted:\n")
cat("Original:\n")
print(eff_size(emm_original, sigma = sigma(best_model), edf = df.residual(best_model)))
cat("\nAge-Restricted:\n")
print(eff_size(emm_restricted, sigma = sigma(best_age_restricted), 
               edf = df.residual(best_age_restricted)))

# Age restirctied plots
p <- ggplot(data_age_restricted, aes(x = eeg_age, y = l_index_frontal, color = affected_group)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=5), 
              aes(fill = affected_group), alpha = 0.2) +
  labs(title = "Lateralization Trajectories by Group",
       x = "Age (years)", y = "Lateralization Index",
       color = "Group", fill = "Group") +
  theme_minimal(base_size = 14)

print(p)

# Save plot
ggsave("output/figures/trajectories_age_restricted.png",  bg = "white",
       width = 10, height = 6, dpi = 300)


p_sex_facet <- ggplot(data_age_restricted, 
                      aes(x = eeg_age, y = l_index_frontal, color = affected_group)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=5), 
              aes(fill = affected_group), alpha = 0.2, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", alpha = 0.6) +

  facet_wrap(~ sex, labeller = labeller(sex = c(female = "Females", male = "Males"))) +
  labs(
    title = "Lateralization Trajectories by Group and Sex",
    x = "Age (years)", 
    y = "Lateralization Index",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

print(p_sex_facet)
ggsave("output/figures/trajectories_by_sex_faceted.png",  bg = "white",
       width = 12, height = 6, dpi = 300)

## SENSITIVITY TWO: FAMILY STRUCTURE


# Identify families with 2+ members
multi_fam_ids <- family_structure %>%
  filter(n_members >= 2) %>%
  pull(fam_id)

data_multi_sibs <- data %>%
  filter(fam_id %in% multi_fam_ids)

# Refit model
multi_sibs <- gam(formula_m3,
  data = data_multi_sibs, REML = TRUE
)

# Compare ICC
icc_original <- performance::icc(best_model)$ICC_adjusted
icc_multi <- performance::icc(multi_sibs)$ICC_adjusted

cat("\nICC Comparison:\n")
cat("Original (with singletons):", round(icc_original, 3), "\n")
cat("Multi-sib families only:", round(icc_multi, 3), "\n")

print(summary(best_model))

# Fit simplified model without random effects
m3_simplified <- gam(l_index_frontal ~ affected_group + 
                       s(eeg_age, by=affected_group, k=5) + sex + site,
                     data = data, method = "REML")

cat("Simplified model (no random effects):\n")
print(summary(m3_simplified))

cat("\n--- AIC COMPARISON ---\n")
cat(sprintf("M3 with family RE: %.2f\n", AIC(m3)))
cat(sprintf("M3 without family RE: %.2f\n", AIC(m3_simplified)))
cat(sprintf("Difference: %.2f\n", AIC(m3_simplified) - AIC(m3)))

# EXPLORATORY GRAPHS 




# Calculate group differences at each age
age_seq <- seq(4, 18, by = 1)
diff_data <- expand.grid(
  eeg_age = age_seq,
  sex = levels(data_model$sex),
  site = "hsj"
)

# Get predictions for each group
for(grp in levels(data_model$affected_group)) {
  diff_data[[grp]] <- predict(m4, 
                              newdata = data.frame(
                                eeg_age = diff_data$eeg_age,
                                affected_group = grp,
                                sex = diff_data$sex,
                                site = diff_data$site
                              ))
}

# Calculate ASD vs Non-Affected difference
diff_data$ASD_vs_NonAff <- diff_data$ASD - diff_data$`Non-Affected`

# Plot
ggplot(diff_data, aes(x = eeg_age, y = ASD_vs_NonAff, 
                      color = sex, linetype = sex)) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(female = "#e91e63", male = "#2196f3"),
                     labels = c("Females", "Males")) +
  scale_linetype_manual(values = c(female = "solid", male = "dashed"),
                        labels = c("Females", "Males")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "ASD Group Difference from Non-Affected Controls",
    subtitle = "Increasing divergence in females, stable difference in males",
    x = "Age (years)",
    y = "Difference in Lateralization Index\n(ASD - Non-Affected)",
    color = "Sex",
    linetype = "Sex"
  ) +
  annotate("text", x = 16, y = max(diff_data$ASD_vs_NonAff) * 0.9,
           label = "Females diverge\nwith age", 
           color = "#e91e63", size = 4, fontface = "bold") +
  annotate("text", x = 16, y = min(diff_data$ASD_vs_NonAff) * 0.5,
           label = "Males show\nstable difference", 
           color = "#2196f3", size = 4, fontface = "bold")


# Check for best model FINAL

# 12.1 Performance metrics

best_model=best_age_restricted
perf <- model_performance(best_model)
print("Overall Model Performance:")
print(perf)

# Save to table
perf_table <- as.data.frame(perf)
write.csv(perf_table, 
          paste0("output/tables/model_performance_", outcome, ".csv"), 
          row.names = FALSE)


# 12.2. R-SQUARED DECOMPOSITION
# Marginal R² = variance explained by FIXED effects only
# Conditional R² = variance explained by FIXED + RANDOM effects

r2_vals <- performance::r2(best_model)
r2_vals

resid_df <- data.frame(
  fitted = fitted(best_model),
  residuals = residuals(best_model),
  eeg_age = data$eeg_age,
  affected_group = data$affected_group,
  sex = data$sex,
  fam_id = data$fam_id
)


# A. Normality of residuals
shapiro_resid <- shapiro.test(resid_df$residuals)
cat("\nShapiro-Wilk test on residuals: W =", round(shapiro_resid$statistic, 4),
    ", p =", round(shapiro_resid$p.value, 4), "\n")

# B. Homoscedasticity
bp_test <- lmtest::bptest(lm(residuals ~ fitted, data = resid_df))
cat("Breusch-Pagan test: χ² =", round(bp_test$statistic, 2),
    ", p =", round(bp_test$p.value, 4), "\n")

# C. Influential observations
check_outliers(best_model)

# D. Multicollinearity (re-check with final model)
check_collinearity(best_model)

# E. Overall diagnostic plot
check_model(best_model, residual_type = "normal")
ggsave("output/diagnostics/model_diagnostics_complete.png",
       width = 12, height = 10, dpi = 300)
 

