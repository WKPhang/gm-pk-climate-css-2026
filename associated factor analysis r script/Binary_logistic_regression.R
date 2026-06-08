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
  "occ_group",
  "education", "past_malaria_bin",
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
    pk_pos = as.factor(pk_pos),
    
    sample_type = as.factor(sample_type),
    site_class = relevel(as.factor(site_class), ref = "Orang Asli"),
    age_group = relevel(as.factor(age_group), ref = "21-30"),
    sex = relevel(as.factor(sex), ref = "Male"),
    race = as.factor(race),
    nationality = as.factor(nationality),
    local = relevel(as.factor(local), ref = "Local"),
    occ_group = as.factor(occ_group),
    education = as.factor(education),
    past_malaria_bin = as.factor(past_malaria_bin),
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
  "education",
  "past_malaria_bin",
  "travel",
  "prevention"
)

#####
# Compare model
#####
glmm_null_1 <- glmmTMB(
  pk_pos ~ 1 + (1 | sampling_site),
  data = my_df,
  family = binomial()
)

glmm_null_2 <- glmmTMB(
  pk_pos ~ 1 + (1 | icemr_full_id),
  data = my_df,
  family = binomial()
)

glmm_null_3 <- glmmTMB(
  pk_pos ~ 1 + (1 | sampling_site/icemr_full_id),
  data = my_df,
  family = binomial()
)

# Check hierachy structure
anova(glmm_null_1, glmm_null_2, glmm_null_3)
icc(glmm_null_1)
icc(glmm_null_2)
icc(glmm_null_3)

glm_null <- glm(
  pk_pos ~ 1,
  data = my_df,
  family = binomial())

glm_f_null <- logistf(
  pk_pos ~ 1,
  data = my_df
)

summary(glm_null)
summary(glm_f_null)

AIC(glmm_null_1, glmm_null_2,glmm_null_3, glm_null, glm_f_null)
AIC(glmm_null_2, glm_null, glm_f_null)
# =========================================================
# 6. Univariable Binary Logistic Regression
# =========================================================
univariable_results <- list()

for (var in predictors) {
  
  fml <- as.formula(paste("pv_pos ~", var))
  
  model <- glm(
    formula = fml,
    data = my_df,
    family = binomial()
  )
  
  tidy_model <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  
  # remove intercept
  tidy_model <- tidy_model %>%
    filter(term != "(Intercept)") %>%
    mutate(predictor = var)
  
  univariable_results[[var]] <- tidy_model
}

# =========================================================
# 7. Combine results
# =========================================================
results_df <- bind_rows(univariable_results) %>%
  select(
    predictor,
    term,
    estimate,      # OR
    conf.low,
    conf.high,
    std.error,
    statistic,
    p.value
  ) %>%
  rename(
    odds_ratio = estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    z_value = statistic
  )

print(results_df)

# =========================================================
# 8. Save results
# =========================================================
write.csv(
  results_df,
  "univariable_binary_logistic_PV_results.csv",
  row.names = FALSE)


# =========================================================
# PK Multivariable Binary Logistic Regression
# =========================================================

multi_model <- glm(
  pk_pos ~ age_group + sex + local,
  data = my_df,
  family = binomial()
)

# =========================================================
# Model summary (log-odds scale)
# =========================================================
summary(multi_model)

vif_values <- car::vif(multi_model)

print(vif_values)


# =========================================================
# Convert to Odds Ratios with 95% CI
# =========================================================
library(broom)

multi_results <- broom::tidy(multi_model, conf.int = TRUE, exponentiate = TRUE)

# Remove intercept for clean table
multi_results <- multi_results[multi_results$term != "(Intercept)", ]

print(multi_results)

# =========================================================
# Save results
# =========================================================
write.csv(
  multi_results,
  "multivariable_logistic_regression_PK_results.csv",
  row.names = FALSE
)


# =========================================================
# PV Multivariable Binary Logistic Regression
# =========================================================

pv_multi_model <- glm(
  pv_pos ~ age_group + sex + local,
  data = my_df,
  family = binomial()
)

# =========================================================
# Model summary (log-odds scale)
# =========================================================
summary(pv_multi_model)

vif_values <- car::vif(pv_multi_model)

print(vif_values)


# =========================================================
# Convert to Odds Ratios with 95% CI
# =========================================================
library(broom)

pv_multi_results <- broom::tidy(pv_multi_model, conf.int = TRUE, exponentiate = TRUE)

# Remove intercept for clean table
pv_multi_results <- pv_multi_results[pv_multi_results$term != "(Intercept)", ]

print(multi_results)

# =========================================================
# Save results
# =========================================================
write.csv(
  pv_multi_results,
  "multivariable_logistic_regression_PV_results.csv",
  row.names = FALSE
)
