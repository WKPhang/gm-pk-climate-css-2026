# Load the readxl package
library(readxl)
library(MASS)
library(ppcor)
library(dplyr)
library(MuMIn)
library(dlnm)
library(ggplot2)
# Set the path to Excel file

file_path <- "PATH/TO/case_data.xlsx"

# Read Sheet 3 and convert to a base R data.frame
full_df <- read_excel(file_path, sheet = "Week_data")

# Check the structure
str(full_df)

colnames(full_df)

# Ensure data is ordered (important for time series)
full_df <- full_df[order(full_df$start_dt), ]

train_df <- full_df[1:574, ]
test_df  <- full_df[575:626, ]

#### Plot ACF and PACF
png(
  filename = "PATH/TO/ACF_PACF_case_counts.png",
  width = 1800,
  height = 900,
  res = 300
)

par(mfrow = c(1, 2))

# ACF
acf(train_df$case_count,
    main = "ACF of Case Counts",
    lag.max = 36)

# PACF
pacf(train_df$case_count,
     main = "PACF of Case Counts",
     lag.max = 36)

par(mfrow = c(1, 1))

dev.off()


###########################################################
# Average Absolute Cross-Correlation Matrix
# Purpose: Lag-based multicollinearity assessment
############################################################

library(reshape2)
library(ggplot2)
library(viridis)
# ----------------------------------------------------------
# 1. Variables of interest (raw data columns)
# ----------------------------------------------------------
vars_ccf <- c("precipitation", "temp_2m", "rh",
              "nino4_smooth", "soi_smooth", "mei_smooth",
              "case_count_l1", "case_count_l2")

data_mat <- train_df[vars_ccf]

# ----------------------------------------------------------
# 2. Function: average absolute cross-correlation
# ----------------------------------------------------------
avg_abs_ccf <- function(x, y, lag_max = 12) {
  
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  ccf_obj <- ccf(x, y, lag.max = lag_max, plot = FALSE)
  
  mean(abs(ccf_obj$acf))
}

# ----------------------------------------------------------
# 3. Initialize matrix
# ----------------------------------------------------------
n <- length(vars_ccf)

ccf_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(ccf_matrix) <- vars_ccf
rownames(ccf_matrix) <- vars_ccf

# ----------------------------------------------------------
# 4. Fill matrix (remove self-correlation)
# ----------------------------------------------------------
for (i in 1:n) {
  for (j in 1:n) {
    
    if (i == j) {
      ccf_matrix[i, j] <- NA
    } else {
      ccf_matrix[i, j] <- avg_abs_ccf(
        data_mat[[i]],
        data_mat[[j]],
        lag_max = 12
      )
    }
    
  }
}

# ----------------------------------------------------------
# 5. Variable labels (for display only)
# ----------------------------------------------------------
var_labels <- c(
  "precipitation" = "Precipitation",
  "temp_2m" = "Temperature",
  "rh" = "Relative Humidity",
  "nino4_smooth" = "Niño 4",
  "soi_smooth" = "SOI",
  "mei_smooth" = "MEI",
  "case_count_l1" = "Cases (1-week lag)",
  "case_count_l2" = "Cases (2-week lag)"
)

# Apply labels AFTER computation
colnames(ccf_matrix) <- var_labels[colnames(ccf_matrix)]
rownames(ccf_matrix) <- var_labels[rownames(ccf_matrix)]

# ----------------------------------------------------------
# 6. Convert to long format
# ----------------------------------------------------------
ccf_long <- melt(ccf_matrix)
colnames(ccf_long) <- c("Var1", "Var2", "AvgAbsCCF")

# ----------------------------------------------------------
# 7. Heatmap (final plot)
# ----------------------------------------------------------
p_ccf_mat<- ggplot(ccf_long, aes(Var1, Var2, fill = AvgAbsCCF)) +
  geom_tile(color = "white") +
  
  # add values on top of tiles
  geom_text(aes(label = round(AvgAbsCCF, 2)), size = 3) +
  
  # viridis color scale
  scale_fill_viridis(option = "C", na.value = "grey90") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  
  labs(
    title = "Average Absolute Cross-Correlation Function Matrix",
    x = NULL,
    y = NULL,
    fill = "Average |CCF|"
  )
p_ccf_mat
ggsave("PATH/TO/p_ccf_mat.png",
       plot = p_ccf_mat,
       width = 8,
       height = 6,
       dpi = 300)


# Null model
null_model <- glm.nb(
  case_count ~ 1 + offset(log(pop_smooth)),
  link = "log",
  data =train_df[13:574,]
)

chat <- sum(residuals(null_model, type = "pearson")^2) /
  df.residual(null_model) #for chat in QAIC,QBIC calculation
chat

#### Crossbasis building #####
# cb.temp (all=2,2; intro_hlc_n=2,3)
# cb.rh (all =2,2)
# cb.prec (all =2,2; second_best = 4,2)

cb.prec <- crossbasis(
  train_df$precipitation,
  lag = 12,
  argvar = list(fun = "ns", df = 5),
  arglag = list(fun = "ns", df = 2)
)

cb.temp <- crossbasis(
  train_df$temp_2m,
  lag = 12,
  argvar = list(fun = "ns", df = 5),
  arglag = list(fun = "ns", df = 2)
)

cb.rh <- crossbasis(
  train_df$rh,
  lag = 12,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 2)
)

cb.nino <- crossbasis(
  train_df$nino4_smooth,
  lag = 12,
  argvar = list(fun = "ns", df = 3),
  arglag = list(fun = "ns", df = 3)
)

cb.soi <- crossbasis(
  train_df$soi_smooth,
  lag = 12,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 3)
)

cb.mei <- crossbasis(
  train_df$mei_smooth,
  lag = 12,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 3)
)


# Crossbasis parameter selection
selected_var <- train_df$mjo_romi
{ # Run on this line -->>>
  
  df_var_options <- 2:5
  df_lag_options <- 2:5
  
  results <- list()
  i <- 1
  
  for (df_var in df_var_options) {
    for (df_lag in df_lag_options) {
      
      cb.selected <- crossbasis(
        selected_var,
        lag = 12,
        argvar = list(fun = "ns", df = df_var),
        arglag = list(fun = "ns", df = df_lag)
      )
      
      formula_selected <- as.formula(
        "case_count ~ cb.selected + offset(log(pop_smooth))"
      )
      
      model <- glm.nb(
        formula_selected,
        link = "log",
        data = train_df
      )
      
      # -----------------------------
      # CONDITION INDEX CALCULATION
      # -----------------------------
      X <- model.matrix(model)
      
      # remove intercept
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      
      # scale for stability (important)
      X_scaled <- scale(X)
      
      eig_vals <- eigen(cor(X_scaled))$values
      
      CN <- sqrt(max(eig_vals) / min(eig_vals))
      
      # extract covariance matrix
      vc <- vcov(model)
      
      # remove intercept if present
      vc <- vc[colnames(vc) != "(Intercept)", colnames(vc) != "(Intercept)"]
      
      # standard errors
      se_vals <- sqrt(diag(vc))
      
      mean_se <- mean(se_vals, na.rm = TRUE)
      max_se  <- max(se_vals, na.rm = TRUE)
      
      
      results[[i]] <- data.frame(
        df_var = df_var,
        df_lag = df_lag,
        AIC = AIC(model),
        AICc = AICc(model),
        BIC = BIC(model),
        condition_number = CN, #>30 = multicollinear
        mean_se = mean_se,
        max_se = max_se
      )
      
      i <- i + 1
    }
  }
  
  results_df <- do.call(rbind, results)
  results_df$diff_AIC <- results_df$AIC - min(results_df$AIC)
  results_df[order(results_df$diff_AIC), ]
}






y <- predict(model) # Sanity check
length(y)
length(train_df$case_count)

# Full model
full_model <- glm.nb(
  case_count ~ cb.prec + cb.temp + cb.rh + 
    cb.nino + cb.soi + cb.mei +
    case_count_l1 + case_count_l2 +  offset(log(pop_smooth)),
  link = "log",
  data = train_df
)

##################
### AIC
#############

# List all probable variables
vars <- c("cb.prec", "cb.temp", "cb.rh", 
          "cb.nino", "cb.soi", "cb.mei",
          "case_count_l1", "case_count_l2")

results <- lapply(vars, function(v) {
  # Reduced model (drop variable)
  reduced <- update(full_model, paste(". ~ . -", v))
  
  # Single-variable model (only that variable)
  single_formula <- as.formula(
    paste("case_count ~", v, "+ offset(log(pop_smooth))")
  )
  
  single <- glm.nb(single_formula, data = train_df)
  
  data.frame(
    variable = v,
    AIC_full = AIC(full_model),
    AIC_reduced = AIC(reduced),
    AIC_single = AIC(single),
    delta_AIC_reduced = AIC(reduced) - AIC(full_model),
    delta_AIC_single = AIC(single) - AIC(full_model) 
  )
})

results_df <- do.call(rbind, results)

# Rank variables (MOST IMPORTANT at bottom, ggplot will revert the order later)
results_df <- results_df[order(results_df$delta_AIC_reduced), ]

label_map <- c(
  "cb.prec" = "Precipitation",
  "cb.temp" = "Temperature",
  "cb.rh"   = "Relative Humidity",
  "cb.nino" = "Niño 4",
  "cb.soi"  = "SOI",
  "cb.mei"  = "MEI",
  "case_count_l1" = "Cases (1-week lag)",
  "case_count_l2" = "Cases (2-week lag)"
)

results_df$variable_label <- label_map[as.character(results_df$variable)]

# 3. Fix factor AFTER sorting (critical step)
results_df$variable_label <- factor(
  results_df$variable_label,
  levels = unique(results_df$variable_label)
)

# subtitle = "Reduced model (leave-one-out) (orange) vs single-variable model (steel blue)"
p_aic <- ggplot(results_df, aes(y = variable_label)) +
  # Back bar (reduced model)
  # Front bar (single-variable)
  geom_segment(aes(x = 0, xend = delta_AIC_reduced),
               yend = results_df$variable_label,
               size = 15, color = "#4e79a7") +
  
  
  xlim(-12, 8) +
  labs(
    x = "ΔAIC",
    y = NULL,
    
    title = "Relative importance based on ΔAIC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11)
  )
p_aic
ggsave("PATH/TO/DLNM_AIC_var_imp.png",
       plot = p_aic,
       width = 8,
       height = 5,
       dpi = 300)

#############
# Deviance loss
#############

# Compute Deviance based importance
results <- lapply(vars, function(v) {
  
  # Reduced model (leave-one-out)
  reduced <- update(full_model, paste(". ~ . -", v))
  
  # Single-variable model
  single_formula <- as.formula(
    paste("case_count ~", v,
          "+ offset(log(pop_smooth))")
  )
  
  single <- glm.nb(single_formula, data = train_df)
  
  data.frame(
    variable = v,
    
    dev_full = deviance(full_model),
    dev_reduced = deviance(reduced),
    dev_single = deviance(single),
    p = anova(full_model, reduced)$`Pr(Chi)`[2], # Log-likelihood test
    # IMPORTANT: deviance LOSS (correct metric)
    delta_dev_reduced = deviance(reduced) - deviance(full_model),
    delta_dev_single  = deviance(single)  - deviance(full_model)
  )
})

results_dev_df <- do.call(rbind, results)

# Rank by importance (largest loss first)
results_dev_df <- results_dev_df[order(- results_dev_df$delta_dev_reduced, decreasing = TRUE), ]

label_map <- c(
  "cb.prec" = "Precipitation",
  "cb.temp" = "Temperature",
  "cb.rh"   = "Relative Humidity",
  "cb.nino" = "Niño 4",
  "cb.soi"  = "SOI",
  "cb.mei"  = "MEI",
  "case_count_l1" = "Cases (1-week lag)",
  "case_count_l2" = "Cases (2-week lag)"
)

results_dev_df$variable_label <- label_map[results_dev_df$variable]

results_dev_df$variable_label <- factor(
  results_dev_df$variable_label,
  levels = results_dev_df$variable_label
)

p_dev_loss <- ggplot(results_dev_df, aes(y = variable_label)) +
  # Reduced model (orange = system loss)
  geom_segment(aes(x = 0, xend = delta_dev_reduced),
               size = 15, color = "#4e79a7") +
  
  # subtitle = "Orange: leave-one-out loss | Blue: single-variable contribution"  
  labs(
    x = "Deviance loss",
    y = NULL,
    title = "Relative importance based on deviance loss",
  ) +
  
  xlim(-4, 3) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11)
  )

p_dev_loss
ggsave("PATH/TO/Dev_loss_var_imp.png",
       plot = p_dev_loss,
       width = 8,
       height = 5,
       dpi = 300)


# Try removing the anomalies and repeat it

##################
### AIC
#############

full_model_2<- glm.nb(
  case_count ~ cb.prec + cb.temp + cb.rh + 
    cb.soi + cb.mjo +
    case_count_l1 + case_count_l2 + offset(log(pop_smooth)),
  link = "log",
  data = train_df
)


# List all probable variables
vars <- c("cb.prec", "cb.temp", "cb.rh",  "cb.soi",
          "case_count_l1", "case_count_l2")

results <- lapply(vars, function(v) {
  # Reduced model (drop variable)
  reduced <- update(full_model, paste(". ~ . -", v))
  
  # Single-variable model (only that variable)
  single_formula <- as.formula(
    paste("case_count ~", v, "+ offset(log(pop_smooth))")
  )
  
  single <- glm.nb(single_formula, data = train_df)
  
  data.frame(
    variable = v,
    AIC_full = AIC(full_model_2),
    AIC_reduced = AIC(reduced),
    AIC_single = AIC(single),
    delta_AIC_reduced = AIC(reduced) - AIC(full_model_2),
    delta_AIC_single = AIC(single) - AIC(full_model_2) 
  )
})

results_df <- do.call(rbind, results)

# Rank variables (MOST IMPORTANT at bottom, ggplot will revert the order later)
results_df <- results_df[order(results_df$delta_AIC_reduced), ]

label_map <- c(
  "cb.prec" = "Precipitation",
  "cb.temp" = "Temperature",
  "cb.rh"   = "Relative Humidity",
  "cb.soi"  = "SOI",
  "case_count_l1" = "Cases (1-week lag)",
  "case_count_l2" = "Cases (2-week lag)"
)

results_df$variable_label <- label_map[as.character(results_df$variable)]

# 3. Fix factor AFTER sorting (critical step)
results_df$variable_label <- factor(
  results_df$variable_label,
  levels = unique(results_df$variable_label)
)

# subtitle = "Reduced model (leave-one-out) (orange) vs single-variable model (steel blue)"
p_aic_2 <- ggplot(results_df, aes(y = variable_label)) +
  # Back bar (reduced model)
  # Front bar (single-variable)
  geom_segment(aes(x = 0, xend = delta_AIC_reduced),
               yend = results_df$variable_label,
               size = 15, color = "#4e79a7") +
  
  xlim(0, 20) +
  labs(
    x = "ΔAIC",
    y = NULL,
    
    title = "Relative importance based on ΔAIC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11)
  )
p_aic_2
ggsave("PATH/TO/DLNM_AIC_2_var_imp.png",
       plot = p_aic_2,
       width = 8,
       height = 5,
       dpi = 300)

#####################
# Final Model Building
#################
# ----------------------------------------------------------
# 1. Fit models (clean indexing via data argument)
# ----------------------------------------------------------
base_model <- glm.nb(
  case_count ~ cb.prec + cb.rh + cb.soi +
    offset(log(pop_smooth)),
  data = train_df,
  link = "log"
)


manipulated_model <- glm.nb(
  case_count ~ cb.prec + cb.rh + cb.soi + case_count_l1 +
    offset(log(pop_smooth)),
  data = train_df,
  link = "log"
)

# ----------------------------------------------------------
# 2. Function to compute evaluation metrics
# ----------------------------------------------------------
model_metrics <- function(model) {
  
  # Predictions
  pred <- predict(model, type = "response")
  obs  <- model$y
  
  # Residuals
  resid <- obs - pred
  
  # RMSE
  rmse <- sqrt(mean(resid^2))
  
  # Predictive R2 (1 - SSE/SST)
  r2 <- 1 - sum(resid^2) / sum((obs - mean(obs))^2)
  
  # AIC, BIC, Deviance
  aic <- AIC(model)
  aicc <- AICc(model)
  bic <- BIC(model)
  dev <- deviance(model)
  
  
  return(c(
    AIC = aic,
    BIC = bic,
    AICc = aicc,
    Deviance = dev,
    RMSE = rmse,
    R2 = r2
  ))
}

# ----------------------------------------------------------
# 3. Compute metrics
# ----------------------------------------------------------
base_metrics <- model_metrics(base_model)
manipulated_metrics <- model_metrics(manipulated_model)

# ----------------------------------------------------------
# 4. Combine into comparison table
# ----------------------------------------------------------
comparison <- rbind(
  Base_Model = base_metrics,
  Manipulated_Model = manipulated_metrics
  
)

comparison
anova(manipulated_model, base_model)$`Pr(Chi)`[2] # Log-likelihood test

manipulated_model <- glm.nb(
  case_count ~ cb.rh + cb.prec + cb.soi + cb.mjo +
    offset(log(pop_smooth)),
  data = train_df,
  link = "log"
)



# Residual diagnostics
res_raw <- residuals(manipulated_model, type = "response")
res_pearson <- residuals(manipulated_model, type = "pearson")
res_dev <- residuals(manipulated_model, type = "deviance")
fitted_vals <- fitted(manipulated_model)


plot(fitted_vals, res_pearson,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, col = "red")

idx <- as.numeric(rownames(model$model))


x <- acf(res_pearson)

plot(fitted(manipulated_model), res_pearson)
abline(h = 0, col = "red")

acf(res_pearson, lag.max = 100)
Box.test(res_pearson, lag = 20, type = "Ljung")
sum(res_pearson^2) / df.residual(manipulated_model)
###### Interpretation contour plot




pred.prec <- crosspred(
  cb.prec, 
  model = manipulated_model,
  by = 0.1,
  cen = mean(full_df$precipitation, na.rm = TRUE),
  model.link = "log"
)


pred.rh <- crosspred(
  cb.rh,
  model = manipulated_model,
  by = 0.1,
  cen = mean(full_df$rh, na.rm = TRUE),
  model.link = "log"
)


pred.soi <- crosspred(
  cb.soi,
  manipulated_model,
  by = 0.01,
  cen = mean(full_df$soi_smooth, na.rm = TRUE),
  model.link = "log"
)


pred.mjo <- crosspred(
  cb.mjo,
  manipulated_model,
  by = 0.01,
  cen = mean(full_df$mjo_romi, na.rm = TRUE),
  model.link = "log"
)

# Help function to snap quantile to nearest grid in crosspred
snap_to_grid <- function(x, grid) {
  grid[which.min(abs(grid - x))]
}


prec_quants <- quantile(full_df$precipitation, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)
prec_quants
prec_quants_matched <- sapply(prec_quants, snap_to_grid, pred.prec$predvar)
prec_quants_matched

rh_quants <- quantile(full_df$rh, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)
rh_quants
rh_quants_matched <- sapply(rh_quants, snap_to_grid, pred.rh$predvar)
rh_quants_matched

soi_quants <- quantile(full_df$soi_smooth, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)
soi_quants
soi_quants_matched <- sapply(soi_quants, snap_to_grid, pred.soi$predvar)
soi_quants_matched

# Proceed to slice plot code for slice plots 

png("PATH/TO/contour_rh.png",
    width = 7, height = 6, units = "in", res = 300)
plot(pred.rh,
     "contour",
     xlab = "Relative humidity (%)",
     ylab = "Lag (weeks)",
     main = "Relative Humidity Exposure–Lag–Response Surface",
     cex.lab = 1.6,   # axis labels
     cex.axis = 1.6,  # tick labels
)
dev.off()

png("PATH/TO/contour_prec.png",
    width = 7, height = 6, units = "in", res = 300)
plot(pred.prec,
     "contour",
     xlab = "24-hr precipitation (mm)",
     ylab = "Lag (weeks)",
     main = "Precipitation Exposure–Lag–Response Surface",
     cex.lab = 1.6,   # axis labels
     cex.axis = 1.6,  # tick labels
)
dev.off()

png("PATH/TO/contour_soi.png",
    width = 7, height = 6, units = "in", res = 300)
contour_soi <- plot(pred.soi,
                    "contour",
                    xlab = "SOI",
                    ylab = "Lag (weeks)",
                    main = "SOI Exposure–Lag–Response Surface",
                    cex.lab = 1.6,   # axis labels
                    cex.axis = 1.6,  # tick labels
)
dev.off()

contour_mjo <- plot(pred.mjo,
                    "contour",
                    xlab = "MJO ROMI",
                    ylab = "Lag (weeks)",
                    main = "MJO ROMI Exposure–Lag–Response Surface",
                    cex.lab = 1.6,   # axis labels
                    cex.axis = 1.6,  # tick labels
)

####################
# FUTURE PREDICTION#
####################

train_df <- full_df[1:574, ]
test_df  <- full_df[575:626, ]

cb.prec_full <- crossbasis(
  full_df$precipitation,
  lag = 12,
  argvar = list(fun = "ns", df = 5),
  arglag = list(fun = "ns", df = 2)
)

cb.rh_full <- crossbasis(
  full_df$rh,
  lag = 12,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 2)
)

cb.soi_full <- crossbasis(
  full_df$soi_smooth,
  lag = 12,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 3)
)

colnames(full_df)

#Rename each block properly
colnames(cb.prec_full) <- paste0("prec_", colnames(cb.prec_full))
colnames(cb.rh_full)   <- paste0("rh_", colnames(cb.rh_full))
colnames(cb.soi_full)  <- paste0("soi_", colnames(cb.soi_full))

merged_df <- cbind(full_df[,c(2,4,12)],
                   as.data.frame(cb.prec_full),
                   as.data.frame(cb.rh_full),
                   as.data.frame(cb.soi_full))

colnames(merged_df)

#Model fitting
train_idx <- 1:574
test_idx  <- 575:626



predictors <- setdiff(colnames(merged_df),
                      c("start_dt","case_count", "pop_smooth"))
formula <- as.formula(
  paste("case_count ~",
        paste(predictors, collapse = " + "),
        "+ offset(log(pop_smooth))")
)

mod_lag0 <- glm.nb(formula, data = merged_df[train_idx, ])

train_count <- predict(
  mod_lag0,
  newdata = merged_df[train_idx, ],
  type = "response"
)

pred_count <- predict(
  mod_lag0,
  newdata = merged_df[test_idx, ],
  type = "response"
)

# Compute BOTH 80% and 95% prediction intervals
# Link-scale prediction
pred_link <- predict(
  mod_lag0,
  newdata = merged_df[test_idx, ],
  type = "link",
  se.fit = TRUE
)

fit_link <- pred_link$fit
se_link  <- pred_link$se.fit

mu <- exp(fit_link)
theta <- mod_lag0$theta

# Variance components
var_mean <- (mu^2) * (se_link^2)
var_nb   <- mu + (mu^2 / theta)
var_total <- var_mean + var_nb
se_total  <- sqrt(var_total)

# Z values
z95 <- 1.96
z80 <- 1.28

# 95% PI
upper_95 <- mu + z95 * se_total
lower_95 <- pmax(0, mu - z95 * se_total)

# 80% PI
upper_80 <- mu + z80 * se_total
lower_80 <- pmax(0, mu - z80 * se_total)

# Create plotting dataframe
plot_df <- merged_df %>%
  mutate(
    observed = case_count,
    fitted = NA,
    predicted = NA,
    lower_95 = NA,
    upper_95 = NA,
    lower_80 = NA,
    upper_80 = NA
  )

# Fill values
plot_df$fitted[train_idx]    <- train_count
plot_df$predicted[test_idx]  <- mu

plot_df$lower_95[test_idx] <- lower_95
plot_df$upper_95[test_idx] <- upper_95
plot_df$lower_80[test_idx] <- lower_80
plot_df$upper_80[test_idx] <- upper_80


# Actual plot
full_ts_plot <- ggplot(plot_df, aes(x = start_dt)) +
  
  # 95% PI (lightest)
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "violet", alpha = 0.3) +
  
  # 80% PI (darker)
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
              fill = "red", alpha = 0.3) +
  
  # observed
  geom_line(aes(y = observed, color = "Observed"), linewidth = 0.5) +
  
  # fitted (train)
  geom_line(aes(y = fitted, color = "Fitted (Train)"), linewidth = 1) +
  
  # predicted (test)
  geom_line(aes(y = predicted, color = "Predicted (Test)"), linewidth = 1) +
  
  # vertical split line
  geom_vline(xintercept = plot_df$start_dt[min(test_idx)],
             linetype = "dashed", color = "gray40") +
  
  scale_color_manual(
    values = c(
      "Observed" = "grey30",
      "Fitted (Train)" = "blue",
      "Predicted (Test)" = "red"
    )
  ) +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    x = "Time",
    y = "Case Count",
    color = ""
  ) +
  
  ylim(0,15) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )
full_ts_plot

ggsave("PATH/TO/Dev_loss_var_imp.png",
       plot = full_ts_plot,
       width = 10,
       height = 5,
       dpi = 300)


# Model performance calculation
# Residual_calculation
resid_train <- merged_df$case_count[train_idx] - train_count
resid_test <- merged_df$case_count[test_idx] - pred_count



# RMSE
rmse_train <- sqrt(mean(resid_train^2, na.rm = T))
rmse_test <- sqrt(mean(resid_test^2, , na.rm = T))


# Final Residual Diagnostic
plot(merged_df$case_count[1:574], train_count,
     xlab = "Observed", ylab = "Predicted")
abline(0, 1, col = "red")

plot(log1p(merged_df$case_count[1:574]), log1p(train_count),
     xlab = "log(Observed + 1)", ylab = "log(Predicted + 1)")
abline(0, 1, col = "red")
