# =========================================================
# 0. Load packages
# =========================================================
library(readxl)
library(dplyr)
library(glmmTMB)
library(performance)
library(boot)
library(MuMIn)
library(performance)
# =========================================================
# 1. Read data
# =========================================================
file_path <- "C:/Users/WeiKitPhang/Documents/UM/Writing/Paper 21 (Gua Musang Preliminary)/Data_processing/factor_analysis.xlsx"

initial_df <- read_excel(file_path, sheet = "analysis_datasheet")
colnames(initial_df)
# =========================================================
# 2. Select variables
# =========================================================
vars <- c(
  "icemr_full_id", "sampling_site", "site_class", "pk_pos", "pv_pos",
  "sample_type", "age_group", "sex", "nationality", "local", "race",
  "occ_group",
  "education", "past_malaria_bin", "fever",
  "travel", "prevention"
)

my_df <- initial_df %>%
  select(all_of(vars)) %>%
  filter(if_all(all_of(vars), ~ !is.na(.)))
nrow(my_df)
# =========================================================
# 3. Convert data types
# =========================================================
my_df <- my_df %>%
  mutate(
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
    fever = as.factor(fever),
    travel = as.factor(travel),
    prevention = as.factor(prevention),
    
    sampling_site = as.factor(sampling_site),
    icemr_full_id = as.factor(icemr_full_id)
  )



# Check if hierachical mixed model is needed
glm_null <- glm(
  pk_pos ~ 1,
  data = my_df,
  family = binomial())


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
icc(glmm_null_1)
icc(glmm_null_2)
icc(glmm_null_3)


performance::check_singularity(glmm_null_1)
performance::check_singularity(glmm_null_2)
performance::check_singularity(glmm_null_3)

# Check for stability
VarCorr(glmm_null_1)

glmmTMB::diagnose(glmm_null_2)

# Singularity check
summary(glmm_null_3)


table_by_id <- my_df %>%
  dplyr::group_by(icemr_full_id) %>%
  dplyr::summarise(
    n = dplyr::n(),
    cases = sum(pk_pos),
    unique_outcome = dplyr::n_distinct(pk_pos)
  )
table(table_by_id$unique_outcome)

table_by_id <- my_df %>%
  dplyr::group_by(icemr_full_id) %>%
  dplyr::summarise(
    n = dplyr::n(),
    cases = sum(pv_pos),
    unique_outcome = dplyr::n_distinct(pv_pos)
  )
table(table_by_id$unique_outcome)




table_by_id <- my_df %>%
  dplyr::group_by(sampling_site) %>%
  dplyr::summarise(
    n = dplyr::n(),
    cases = sum(pk_pos),
    unique_outcome = dplyr::n_distinct(pk_pos)
  )
table(table_by_id$unique_outcome)



# =========================================================
# 4. Final predictor list 
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
  "fever",
  "travel",
  "prevention"
)

# =========================================================
# 5. Univariable GLMMTMB (PK outcome)
# =========================================================

univariable_results <- list()

for (var in predictors) {
  
  # force clean formula construction
  fml <- as.formula(paste("pk_pos ~", var, "+ (1 | icemr_full_id)"))
  
  model <- glmmTMB(
    formula = fml,
    data = my_df,
    family = binomial()
  )
  
  coef_table <- summary(model)$coefficients$cond
  
  # skip intercept row
  coef_table <- coef_table[rownames(coef_table) != "(Intercept)", , drop = FALSE]
  
  OR <- exp(coef_table[, "Estimate"])
  lower <- exp(coef_table[, "Estimate"] - 1.96 * coef_table[, "Std. Error"])
  upper <- exp(coef_table[, "Estimate"] + 1.96 * coef_table[, "Std. Error"])
  
  results <- data.frame(
    predictor = var,
    term = rownames(coef_table),
    odds_ratio = OR,
    conf_low = lower,
    conf_high = upper,
    std_error = coef_table[, "Std. Error"],
    z_value = coef_table[, "z value"],
    p_value = coef_table[, "Pr(>|z|)"]
  )
  
  univariable_results[[var]] <- results
}

# =========================================================
# 6. Combine results
# =========================================================
results_df <- bind_rows(univariable_results)

print(results_df)

# =========================================================
# 7. Save results
# =========================================================
write.csv(results_df,
          "univariable_GLMMTMB_PK_results.csv",
          row.names = FALSE)


# =========================================================
# PK Multivariable Binary Logistic Regression
# =========================================================

multi_model <- glmmTMB(
  pk_pos ~ age_group + sex + fever + (1 | icemr_full_id),
  data = my_df,
  family = binomial()
)

summary(multi_model)
VarCorr(multi_model)
