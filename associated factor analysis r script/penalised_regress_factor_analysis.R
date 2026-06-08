# =========================================================
# 0. Load packages
# =========================================================
library(readxl)
library(dplyr)
library(broom)
library(purrr)
library(car)
library(glmmTMB)
library(MuMIn)
library(logistf)
library(detectseparation)

?vcovCL

# =========================================================
# 1. Read data
# =========================================================
file_path <- "C:/Users/WeiKitPhang/Documents/UM/Writing/Paper 21 (Gua Musang Preliminary)/Data_processing/factor_analysis.xlsx"

initial_df <- read_excel(file_path, sheet = "analysis_datasheet")

# =========================================================
# 2. Select variables
# =========================================================
vars <- c(
  "icemr_full_id", "sampling_site", "site_class", "pk_pos", "pv_pos",
  "sample_type", "age_group", "sex", "nationality", "local", "race",
  "occ_group", "agri_work",
  "education", "past_malaria_bin", "fever",
  "travel", "prevention"
)

my_df <- initial_df %>%
  select(all_of(vars)) %>%
  filter(if_all(all_of(vars), ~ !is.na(.)))

# =========================================================
# 3. Convert data types
# =========================================================
my_df <- my_df %>%
  mutate(
    sample_type = as.factor(sample_type),
    sampling_site = relevel(as.factor(sampling_site), ref = "Panggong Lalat"),
    site_class = relevel(as.factor(site_class), ref = "Orang Asli"),
    age_group = relevel(as.factor(age_group), ref = "21-30"),
    sex = relevel(as.factor(sex), ref = "Male"),
    race = relevel(as.factor(race), ref = "Orang Asli"),
    nationality = relevel(as.factor(nationality), ref = "Malaysia"),
    local = relevel(as.factor(local), ref = "Local"),
    occ_group = as.factor(occ_group),
    agri_work = relevel(as.factor(agri_work), ref = "non_agri"),
    education = as.factor(education),
    past_malaria_bin = as.factor(past_malaria_bin),
    fever = as.factor(fever),
    travel = as.factor(travel),
    prevention = as.factor(prevention),
    
    sampling_site = as.factor(sampling_site),
    icemr_full_id = as.factor(icemr_full_id)
  )

# =========================================================
# 5. Predictor list
# =========================================================
predictors <- c(
  "sample_type",
  "sampling_site",
  "site_class",
  "age_group",
  "sex",
  "race",
  "nationality",
  "local",
  "occ_group",
  "agri_work",
  "education",
  "past_malaria_bin",
  "fever",
  "travel",
  "prevention"
)

# =========================================================
# 6. Univariable Binary Logistic Regression
# =========================================================
univariable_results <- list()

for (var in predictors) {
  
  fml <- as.formula(paste("pv_pos ~", var))
  
  model <- logistf(
    formula = fml,
    data = my_df
  )
  
  # Extract results manually
  res <- data.frame(
    term      = names(coef(model)),
    estimate  = coef(model),
    std.error = sqrt(diag(model$var)),
    p.value   = model$prob
  )
  
  # Add CI (log scale → exponentiate later)
  ci <- confint(model)
  res$conf.low  <- ci[, 1]
  res$conf.high <- ci[, 2]
  
  # Convert to odds ratio
  res <- res %>%
    mutate(
      estimate  = exp(estimate),
      conf.low  = exp(conf.low),
      conf.high = exp(conf.high),
      predictor = var
    ) %>%
    filter(term != "(Intercept)")
  
  univariable_results[[var]] <- res
}

# Combine all results
final_results <- bind_rows(univariable_results)


# =========================================================
# 8. Save results
# =========================================================
write.csv(
  final_results,
  "univariable_penalized_logistic_PV_results.csv",
  row.names = FALSE)


# =========================================================
# PK Multivariable Binary Logistic Regression
# =========================================================

multi_model <- logistf(
 pk_pos ~ age_group + sex +
    site_class,
  data = my_df
)

# ---- Extract results correctly ----
results <- data.frame(
  variable = names(multi_model$coefficients),
  estimate = multi_model$coefficients,
  lower_CI = multi_model$ci.lower,
  upper_CI = multi_model$ci.upper,
  p_value = multi_model$prob
)

# ---- Compute OR and CI ----
results <- results %>%
  mutate(
    OR = exp(estimate),
    OR_lower = exp(lower_CI),
    OR_upper = exp(upper_CI),
    significant = p_value < 0.05
  )

# ---- Create clean OR (95% CI) column ----
results <- results %>%
  mutate(
    OR_CI = sprintf("%.2f (%.2f–%.2f)", OR, OR_lower, OR_upper)
  )

# ---- Optional: remove intercept ----
results_clean_Pk <- results %>%
  filter(variable != "(Intercept)")

# ---- View full table ----
results_clean_Pk

# ---- View only significant variables ----
significant_results <- results_clean %>%
  filter(significant == TRUE)

significant_results

# =========================================================
# Detect separation
# =========================================================
glm_model <- glm(
  pk_pos ~ age_group + sex + sampling_site +local,
  data = my_df,
  family = binomial()
)

summary(glm_model)

# =========================================================
# Log-likelihood test
# =========================================================
# IF  p < 0.05, keep the variable
drop1(multi_model, test = "Chisq")


model_reduced <- logistf(
  pk_pos ~ age_group + sex,
  data = my_df
)

# Compare
# < 0.05 full model is better
AIC(multi_model, model_reduced)

anova(model_reduced, multi_model)

#Contingency table
table(my_df$age_group, my_df$pk_pos)
table(my_df$sex, my_df$pk_pos)
table(my_df$site_class, my_df$pk_pos)
# =========================================================
# Save results
# =========================================================
write.csv(
  results_clean_Pk,
  "multivariable_penalized_logistic_PK_results.csv",
  row.names = FALSE
)


# =========================================================
# PV Multivariable Binary Logistic Regression
# =========================================================

multi_model <- logistf(
  pv_pos ~ age_group + sex + local,
  data = my_df
)

# Check VIF using GLM if got issues
glm_temp <- glm(pv_pos ~ age_group + sex + race + agri_work,
                data = my_df, family = binomial())

vif(glm_temp)

# Check contingency table
table(my_df$age_group, my_df$pv_pos)
table(my_df$sex, my_df$pv_pos)
table( my_df$site_class, my_df$pv_pos)
table( my_df$race, my_df$pv_pos)
table( my_df$nationality, my_df$pv_pos)
table( my_df$occ_group, my_df$pv_pos)


# ---- Extract results correctly ----
results <- data.frame(
  variable = names(multi_model$coefficients),
  estimate = multi_model$coefficients,
  lower_CI = multi_model$ci.lower,
  upper_CI = multi_model$ci.upper,
  p_value = multi_model$prob
)

# ---- Compute OR and CI ----
results <- results %>%
  mutate(
    OR = exp(estimate),
    OR_lower = exp(lower_CI),
    OR_upper = exp(upper_CI),
    significant = p_value < 0.05
  )

# ---- Create clean OR (95% CI) column ----
results <- results %>%
  mutate(
    OR_CI = sprintf("%.2f (%.2f–%.2f)", OR, OR_lower, OR_upper)
  )

# ---- Optional: remove intercept ----
results_clean_Pv <- results %>%
  filter(variable != "(Intercept)")

# ---- View full table ----
results_clean_Pv

# ---- View only significant variables ----
significant_results <- results_clean_Pv %>%
  filter(significant == TRUE)

significant_results

# =========================================================
# Detect separation
# =========================================================
glm_model <- glm(
  pk_pos ~ age_group + sex + sampling_site +local,
  data = my_df,
  family = binomial()
)

summary(glm_model)

# =========================================================
# Log-likelihood test
# =========================================================
# IF  p < 0.05, keep the variable
drop1(multi_model, test = "Chisq")


model_reduced <- logistf(
  pv_pos ~ age_group + sex + site_class,
  data = my_df
)

# Compare
# < 0.05 full model is better
AIC(multi_model, model_reduced)

anova(model_reduced, multi_model)

#Contingency table
table(my_df$age_group, my_df$pv_pos)
table(my_df$sex, my_df$pv_pos)
table(my_df$race, my_df$pv_pos)
# =========================================================
# Save results
# =========================================================
write.csv(
  results_clean_Pv,
  "multivariable_penalized_logistic_PV_results.csv",
  row.names = FALSE
)

