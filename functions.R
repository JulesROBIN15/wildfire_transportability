#######################################################################################################
### Data loading
#######################################################################################################

#' Load and prepare wildfire smoke exposure and demographic data
#'
#' @description
#' Downloads wildfire smoke exposure data by California county (2018), 
#' calculates average exposure during the risk window, integrates census 
#' demographic population data, and constructs time-indexed exposure 
#' variables for analyses.
#'
#' @return A data.frame containing county-level daily observations with the 
#'   following variables:
#'   - COUNTY_1 : Numeric county identifier ;
#'   - year, month, day : Date components of the observation ;
#'   - week2 : Formatted week string ("2018-W37") ;
#'   - countyname : Name of the California county ;
#'   - smoke : Daily wildfire smoke exposure level ;
#'   - tmax : Maximum daily temperature ;
#'   - prec : Daily precipitation ;
#'   - hum : Relative humidity ;
#'   - shrtwv_rad : Shortwave radiation ;
#'   - wind : Wind speed ;
#'   - circ : Cardiovascular hospital admissions count ;
#'   - resp : Respiratory hospital admissions count ;
#'   - week : Numeric week index ;
#'   - smoke_week45_49 : Mean smoke exposure during the risk window (weeks 45–49) ;
#'   - Smokebin : Binary indicator of high exposure (1 if >= 10, 0 otherwise) ;
#'   - total.pop : Total county population from ACS estimates ;
#'   - time_unit : Sequential time index per county ;
#'   - exp_day : Binary exposure indicator (Smokebin == 1 & week >= 45) ;
#'   - post_time : Elapsed time relative to exposure onset (time_unit - 56).
#'

read_wildfire_data <- function() {
  
  # Load and prepare wildfire data
  CA_hosp_County <- read.csv("https://raw.githubusercontent.com/benmarhnia-lab/Wildfires_social_vulnerability/refs/heads/main/CA_hosp_County_2018.csv")
  
  # Define study period: from Sept 13 to Dec 5, 2018
  CA_hosp_County_week37_49 <- CA_hosp_County[which(CA_hosp_County$week>36 & CA_hosp_County$week<49),]
  
  # Define exposure period: from Nov 8 to Dec 5, 2018 
  CA_hosp_County_week45_49<- CA_hosp_County[which(CA_hosp_County$week>44 & CA_hosp_County$week<49),]
  
  # Average smoke exposure during exposure period
  Counties_exp <- CA_hosp_County_week45_49 |>
    group_by(COUNTY_1) |>
    summarize(smoke_week45_49=mean(smoke))
  
  CA_hosp_County_week37_49 <- left_join(CA_hosp_County_week37_49, Counties_exp,
                                        by = c("COUNTY_1"))
  
  # Define smoke exposure level: average smoke exposure greater than 10
  CA_hosp_County_week37_49$Smokebin<-ifelse(CA_hosp_County_week37_49$smoke_week45_49>=10, 1,
                                            ifelse(CA_hosp_County_week37_49$smoke_week45_49<10, 0, NA))
  
  CA_hosp_County_week37_49$week<- as.numeric(CA_hosp_County_week37_49$week)
  
  ### Load and prepare demographic data: extract total population by county
  dp5 <- read.csv("ACSDP5Y2020.DP05-2025-02-19T205730.csv")
  dp5 <- dp5[2, colnames(dp5)[grepl("Estimate", colnames(dp5))]]
  colnames(dp5) <- str_extract(colnames(dp5), regex("(^.+)County"))
  colnames(dp5) <- str_replace(colnames(dp5), ".County", "")
  colnames(dp5) <- str_replace_all(colnames(dp5), "\\.", " ")
  dp5 <- data.frame(t(dp5))
  dp5 <- cbind(rownames(dp5), dp5)
  rownames(dp5) <- NULL
  colnames(dp5) <- c("countyname", "total.pop")
  dp5$total.pop <- as.numeric(gsub(",", "", dp5$total.pop))
  
  ### Combine exposure and demographic datasets
  data <- left_join(CA_hosp_County_week37_49, dp5, by = c("countyname"))
  data$time_unit <- rep(1:84, 58)
  data$exp_day <- ifelse(data$Smokebin==1 & data$week >= 45, 1, 0)
  data$post_time <- ifelse(data$exp_day==1, data$time_unit - 56, 0)
  return(data)
}



#######################################################################################################
### Effect estimation
#######################################################################################################

#' Estimate wildfire smoke exposure effect using interrupted time series regression
#'
#' @description
#' Fits a generalized linear model (quasi-Poisson for counts with population offset, 
#' or Gaussian for proportions) on subsetted county data to estimate the effect of 
#' smoke exposure, controlling for underlying temporal trends and meteorological
#' confounders.
#'
#' @param county_ind Logical vector indexing rows corresponding to a specific county.
#' @param count Logical. If TRUE (default), fits a quasi-Poisson regression on
#'   respiratory outcome counts with a log population offset. If FALSE, fits
#'   a Gaussian regression on outcome proportions.
#'
#' @return A numeric vector of length 3 containing the point estimate and its 
#'   confidence interval bounds.

estimate_effect <- function(county_ind, count = TRUE) {
  
  # Interrupted time series regression analysis
  if (count == TRUE) {
    
    county_model <- glm(resp ~ offset(log(total.pop)) + 
                          exp_day + 
                          time_unit +
                          tmax + prec + hum + shrtwv_rad + wind,
                        family=quasipoisson,
                        data=exposed_data[county_ind, ])
  } else {
    
    county_model <- glm(resp_prop ~ exp_day + 
                          time_unit +
                          tmax + prec + hum + shrtwv_rad + wind,
                        family=gaussian,
                        data=exposed_data[county_ind, ])
  }
  
  # Isolate point estimate and confidence interval
  estimate <- county_model$coefficients[2][1]
  estimate_ci <- confint(county_model)[2, ]
  
  if (count == TRUE) {
    estimate <- exp(estimate)
    estimate_ci <- exp(estimate_ci)
  }
  
  return(c(estimate, estimate_ci))
}



#######################################################################################################
### Causal transportability
#######################################################################################################

#' Transport effect estimates between counties using Inverse Odds of Sampling Weights (IOSW)
#'
#' @description
#' Estimates transported causal effects from source counties to target counties 
#' using Inverse Odds of Sampling Weights (IOSW). Sampling propensity scores 
#' (county inclusion probabilities) are estimated via both parametric logistic 
#' regression and SuperLearner.
#'
#' @param county_1 Character vector of target county names to transport the effect from.
#' @param county_0 Character vector of source county names to transport the effect to.
#' @param zsite.SL.library Character vector specifying SuperLearner algorithms. 
#'   Defaults to c('SL.glm', 'SL.ranger', 'SL.glmnet').
#'
#' @return A tibble containing transported effect estimates (RR) and confidence
#'   intervals for both logistic regression and SuperLearner methods.

transport_iosw <- function(county_1, county_0, 
                           zsite.SL.library = c('SL.glm', 'SL.ranger', 'SL.glmnet')) {
  
  # Define indicator for counties to transport to/from
  county_1_ind <- which(exposed_data$countyname %in% county_1)
  county_0_ind <- which(exposed_data$countyname %in% county_0)
  
  # Construct the new dataset
  county_1_data <- exposed_data |>
    mutate(COUNTY_1 = case_when(countyname %in% county_1 ~ 1,
                                countyname %in% county_0 ~ 0,
                                T ~ NA)) |>
    filter(!is.na(COUNTY_1))
  
  # Estimate the probability of county inclusion using logistic regression
  cps_lr <- predict(glm(COUNTY_1 ~ tmax + prec + hum + shrtwv_rad + wind,
                        family = binomial,
                        data = county_1_data),
                    type = "response")
  
  # Estimate the probability of county inclusion using SuperLearner
  cps_sl <- SuperLearner(Y=county_1_data$COUNTY_1,
                         X=county_1_data[, c("tmax", "prec", "hum", "shrtwv_rad", "wind")],
                         SL.library=zsite.SL.library,
                         family='binomial')$SL.predict[,1]
  
  # Calculate Inverse Odds of Sampling Weights (IOSW)
  p_county_1 <- sum(county_1_data$COUNTY_1 == 1)/nrow(county_1_data)
  p_county_0 <- sum(county_1_data$COUNTY_1 == 0)/nrow(county_1_data)
  
  county_1_data <- county_1_data |>
    mutate(iosw_lr = ifelse(county_1_data$COUNTY_1 == 1, (1 - cps_lr)/cps_lr * p_county_1/p_county_0, 0),
           iosw_sl = ifelse(county_1_data$COUNTY_1 == 1, (1 - cps_sl)/cps_sl * p_county_1/p_county_0, 0))
  
  # county_1_data <- county_1_data |>
  #     group_by(COUNTY_1) |>
  #     mutate(iosw = pmin(pmax(iosw, quantile(iosw, 0.05)), quantile(iosw, 0.95)))
  
  # Final outcome models
  outcome_mod_lr <- glm(resp ~ offset(log(total.pop)) + 
                          exp_day + 
                          time_unit + 
                          tmax + prec + hum + shrtwv_rad + wind, 
                        family = quasipoisson, 
                        data = county_1_data, 
                        weights = county_1_data$iosw_lr,
                        subset = county_1_data$iosw_lr != 0)
  
  outcome_mod_sl <- glm(resp ~ offset(log(total.pop)) + 
                          exp_day + 
                          time_unit + 
                          tmax + prec + hum + shrtwv_rad + wind, 
                        family = quasipoisson, 
                        data = county_1_data, 
                        weights = county_1_data$iosw_sl,
                        subset = county_1_data$iosw_sl != 0)
  
  # Isolate transported estimates and confidence intervals
  estimate_lr <- exp(outcome_mod_lr$coefficients[2][1])
  ci_lr <- exp(confint(outcome_mod_lr)[2,])
  estimate_sl <- exp(outcome_mod_sl$coefficients[2][1])
  ci_sl <- exp(confint(outcome_mod_sl)[2,])
  
  transport_dat <- tibble(type = c("Transported (logistic regression)", "Transported (SuperLearner)"),
                          estimate = c(estimate_lr, estimate_sl),
                          lower_ci = c(ci_lr[1], ci_sl[1]),
                          upper_ci = c(ci_lr[2], ci_sl[2]))
  
  return(transport_dat)
}



#######################################################################################################
### Visualizing transported effects
#######################################################################################################

#' Plot transported effect estimates from individual and joint source counties to a target county
#'
#' @description
#' Generates a ggplot2 visualization comparing transported effect estimates (rate ratios) 
#' and their 95% confidence intervals across individual source counties (pairwise) 
#' and jointly from all source counties combined to a target county, stratified by 
#' weighting method (e.g., SuperLearner vs. GLM).
#'
#' @param transport_pairwise A data.frame or tibble containing pairwise transported
#'    effect estimates, confidence intervals, source county names, and weighting types.
#' @param transport_jointly A data.frame or tibble containing joint transported
#'    effect estimates (combining all source counties), confidence intervals,
#'    and weighting types.
#' @param county_0 Character string specifying the target county name.
#' @param title Character string specifying the plot title. Optional (defaults to NULL).
#'
#' @return A ggplot object displaying pairwise and joint transported effect estimates
#'  with their corresponding confidence intervals.

plot_transported_effects <- function(transport_pairwise, transport_jointly, county_0, title = NULL) {
  
  # Combine pairwise and joint transport estimates
  combined_transport <- bind_rows(
    transport_pairwise,
    transport_jointly[2:3, ] |> mutate(county = paste0("All of them to ", county_0))
  ) |> 
    mutate(
      county = factor(county, levels = c(unique(transport_pairwise$county), 
                                         paste0("All of them to ", county_0)))
    )
  
  pd <- position_dodge(width = 0.5)
  custom_colors <- c("#F8766D", "#7570B3", "#E69F00")
  
  graph <- combined_transport |>
    ggplot(aes(x = county, y = estimate, color = type, linetype = type)) +
    geom_vline(
      xintercept = length(unique(transport_pairwise$county)) + 0.5,
      linetype = "dotted",
      color = "grey60",
      linewidth = 0.6
    ) +
    geom_errorbar(
      aes(ymin = lower_ci, ymax = upper_ci),
      linewidth = 0.8,
      width = 0.25,
      position = pd
    ) +
    geom_point(
      size = 1.8,
      position = pd
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(breaks = seq(-10, 10, by = 0.2)) +
    labs(
      x = NULL,
      y = "Effect Estimate (RR)",
      color = "Estimate",
      linetype = "Estimate",
      title = title
    ) + 
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.title.y = element_text(size = 10),
      panel.grid.major.x = element_blank()
    )
  
  return(graph)
}



#' Run transportability analysis and plot source-to-target county estimates
#'
#' @description
#' Computes pairwise and joint transported effect estimates using Inverse Odds of 
#' Sampling Weights (IOSW) from a set of source counties to a target county, 
#' combines them with the target county's observed effect, and displays 
#' the comparison plot.
#'
#' @param county_0 Character string specifying the target county name.
#' @param counties_1 Character vector of source county names to transport effects from.
#' @param title Character string specifying the title for the generated plot. 
#'    Optional (defaults to NULL).

transport_analysis <- function(county_0, counties_1, title = NULL) {
  
  observed <- filter(exposed_county_effects, county == county_0) |>
    mutate(type = "Observed")
  
  # Estimate pairwise transported effects from each source county
  transport_pairwise <- t(sapply(counties_1, function(s) transport_iosw(county_1 = s, county_0 = county_0)))
  
  transport_pairwise <- as.data.frame(apply(transport_pairwise, 2, unlist)) |>
    mutate(county = paste0(rep(rownames(transport_pairwise), each = 2), " to ", county_0),
           across(c(estimate, lower_ci, upper_ci), as.numeric)) |>
    bind_rows(observed)
  
  # Estimate joint transported effect from all source counties combined
  transport_jointly <- transport_iosw(counties_1, county_0)
  
  transport_jointly <- observed |>
    dplyr::select(-county) |>
    bind_rows(transport_jointly)
  
  # Generate and display comparison plot
  print(plot_transported_effects(transport_pairwise, transport_jointly, county_0, title))
}


#######################################################################################################
### Covariate Overlap Tools
#######################################################################################################

#' Calculate Standardized Mean Difference (SMD)
#'
#' @description
#' Computes the unweighted standardized mean difference for a continuous variable
#' between two groups (source vs. target county).
#'
#' @param x Numeric vector of covariate values.
#' @param g Binary vector indicating group membership (1 for source, 0 for target).
#'
#' @return Numeric value representing the unweighted SMD.

calc_smd <- function(x, g) {
  
  m1 <- mean(x[g == 1], na.rm = TRUE)
  m0 <- mean(x[g == 0], na.rm = TRUE)
  
  s1 <- var(x[g == 1], na.rm = TRUE)
  s0 <- var(x[g == 0], na.rm = TRUE)
  
  return((m1 - m0) / sqrt((s1 + s0) / 2))
}



#' Calculate Weighted Standardized Mean Difference (WSMD)
#'
#' @description
#' Computes the weighted standardized mean difference for a continuous variable 
#' between two groups, applying weights (IOSW) to group 1 (source).
#'
#' @param x Numeric vector of covariate values.
#' @param g Binary vector indicating group membership (1 for source, 0 for target).
#' @param w Numeric vector of sampling or transport weights applied to group 1.
#'
#' @return Numeric value representing the weighted SMD.

calc_wsmd <- function(x, g, w) {
  
  m1 <- weighted.mean(x[g == 1], w[g == 1], na.rm = TRUE)
  m0 <- mean(x[g == 0], na.rm = TRUE)
  
  s1 <- Hmisc::wtd.var(x[g == 1], weights = w[g == 1], na.rm = TRUE)
  s0 <- var(x[g == 0], na.rm = TRUE)
  
  return((m1 - m0) / sqrt((s1 + s0) / 2))
}



#' Calculate Tipton's Overlap Index (Beta)
#'
#' @description
#' Computes Tipton's Beta index measuring the degree of overlap between 
#' two propensity score distributions (ranging from 0 to 1). A value closer to 1 
#' indicates strong distribution overlap, facilitating causal transportability.
#'
#' @param p Numeric vector of estimated selection probabilities (propensity scores).
#' @param g Binary vector indicating group membership (1 for source, 0 for target).
#' @param common_support Logical. If TRUE, restricts integration to the region 
#'   of shared common support. Default is TRUE.
#' @param n Integer specifying the number of evaluation points for density estimation. 
#'   Default is 1024.
#'
#' @return Numeric value between 0 and 1 representing Tipton's Beta overlap statistic.

calc_tipton_beta <- function(p, g, common_support = TRUE, n = 1024) {
  
  # Extract the inclusion probabilities for each group
  p1 <- p[g == 1]
  p0 <- p[g == 0]
  
  if (common_support) {
    
    # Study only on the common support 
    a <- max(min(p1, na.rm = TRUE), min(p0, na.rm = TRUE))
    b <- min(max(p1, na.rm = TRUE), max(p0, na.rm = TRUE))
    
    if (a >= b) return(0) # No common support
    
    # Estimate the probability of being in the common support.
    p1_C <- mean(p1 >= a & p1 <= b, na.rm = TRUE)
    p0_C <- mean(p0 >= a & p0 <= b, na.rm = TRUE)
    
  } else {
    
    a <- min(c(p1, p0), na.rm = TRUE)
    b <- max(c(p1, p0), na.rm = TRUE)
    
    p1_C <- p0_C <- 1
  }
  
  # Estimate the propensity score distributions between a and b
  d1 <- density(p1, from = a, to = b, n = n, na.rm = TRUE)
  d0 <- density(p0, from = a, to = b, n = n, na.rm = TRUE)
  
  f1 <- pmax(0, d1$y) / p1_C
  f0 <- pmax(0, d0$y) / p0_C
  
  # Compute Tipton's Overlap Index
  dx <- d1$x[2] - d1$x[1]
  integrand <- sqrt(f1 * f0)
  beta <- sqrt(p1_C * p0_C) * sum((integrand[-1] + integrand[-length(integrand)]) / 2) * dx
  
  return(beta)
}


#' Plot Propensity Score Densities and Tipton's Beta
#'
#' @description
#' Fits Logistic Regression and SuperLearner models to estimate selection scores,
#' calculates Tipton's Beta overlap indices, and generates density plots.
#'
#' @param plot_data A data.frame containing environmental covariates and source county 
#'    indicator.
#' @param zsite.SL.library Character vector of SuperLearner algorithms.
#' @param vars Character vector of covariate column names used for modeling selection.
#'
#' @return A data.frame augmented with propensity score predictions (cps_lr, cps_sl).

plot_ps_overlap <- function(plot_data, zsite.SL.library, vars) {
  
  # Fit selection propensity score via Logistic Regression
  formula_lr <- as.formula(paste("COUNTY_1 ~", paste(vars, collapse = " + ")))
  cps_lr <- predict(glm(formula_lr, family = binomial, data = plot_data), type = "response")
  
  # Fit selection propensity score via SuperLearner
  sl_fit <- SuperLearner(
    Y = plot_data$COUNTY_1,
    X = plot_data[, vars],
    SL.library = zsite.SL.library,
    family = 'binomial'
  )
  cps_sl <- sl_fit$SL.predict[, 1]
  
  plot_data <- plot_data |>
    mutate(
      cps_lr = cps_lr,
      cps_sl = cps_sl,
      County_Group = factor(ifelse(COUNTY_1 == 1, "Source", "Target"))
    )
  
  # Compute Tipton's Beta
  beta_lr_global <- calc_tipton_beta(plot_data$cps_lr, plot_data$COUNTY_1, common_support = FALSE)
  beta_lr_cs     <- calc_tipton_beta(plot_data$cps_lr, plot_data$COUNTY_1, common_support = TRUE)
  beta_sl_global <- calc_tipton_beta(plot_data$cps_sl, plot_data$COUNTY_1, common_support = FALSE)
  beta_sl_cs     <- calc_tipton_beta(plot_data$cps_sl, plot_data$COUNTY_1, common_support = TRUE)
  
  # Build density plots
  p_lr <- ggplot(plot_data, aes(x = cps_lr, fill = County_Group)) +
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("#5B2C6F", "#A569BD")) +
    theme_minimal() +
    labs(
      title = sprintf("PS Densities - Logistic Regression (Beta Global = %.3f | Beta CS = %.3f)", beta_lr_global, beta_lr_cs),
      x = "Probability of selection (propensity score)",
      y = "Density",
      fill = "Counties"
    ) +
    theme(legend.position = "bottom")
  
  p_sl <- ggplot(plot_data, aes(x = cps_sl, fill = County_Group)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("#D58000","#F2B222")) +
    theme_minimal() +
    labs(
      title = sprintf("PS Densities - SuperLearner (Beta Global = %.3f | Beta CS = %.3f)", beta_sl_global, beta_sl_cs),
      x = "Probability of selection (propensity score)",
      y = "Density",
      fill = "Counties"
    ) +
    theme(legend.position = "bottom")
  
  print(p_lr)
  print(p_sl)
  
  return(plot_data)
}


#' Love plot before and after IOSW
#'
#' @description
#' Calculates Inverse Odds of Sampling Weights (IOSW) and visualizes standardized mean 
#' differences (SMD) across environmental covariates before and after weighting.
#'
#' @param plot_data A data.frame output from plot_ps_overlap containing cps_lr and cps_sl.
#' @param vars Character vector of covariate column names to evaluate balance for.
#'  Default is c("tmax", "prec", "hum", "shrtwv_rad", "wind").
#' @param var_labels Named character vector mapping raw variable names to new labels. 
#'  Default is c("tmax" = "Max. Temperature", "prec" = "Precipitation", 
#'  "hum" = "Relative Humidity", "shrtwv_rad" = "Shortwave Radiation",
#'  "wind" = "Wind Speed")
#'  
#' @return Displays a love plot showing SMD reduction across covariates.


plot_covariate_balance <- function(plot_data, 
                                   vars = c("tmax", "prec", "hum", "shrtwv_rad", "wind"),
                                   var_labels = c(
                                     "tmax"       = "Max. Temperature",
                                     "prec"       = "Precipitation",
                                     "hum"        = "Relative Humidity",
                                     "shrtwv_rad" = "Shortwave Radiation",
                                     "wind"       = "Wind Speed"
                                   )) {
  
  p_county_1 <- sum(plot_data$COUNTY_1 == 1) / nrow(plot_data)
  p_county_0 <- sum(plot_data$COUNTY_1 == 0) / nrow(plot_data)
  
  # Compute IOSW weights
  plot_data <- plot_data |>
    mutate(
      iosw_lr = ifelse(COUNTY_1 == 1, (1 - cps_lr) / cps_lr * p_county_1 / p_county_0, 0),
      iosw_sl = ifelse(COUNTY_1 == 1, (1 - cps_sl) / cps_sl * p_county_1 / p_county_0, 0)
    )
  
  # Compute SMDs
  smd_before <- data.frame(
    variable = vars,
    smd = sapply(vars, function(v) abs(calc_smd(plot_data[[v]], plot_data$COUNTY_1))),
    stage = "Unweighted"
  )
  
  smd_lr <- data.frame(
    variable = vars,
    smd = sapply(vars, function(v) abs(calc_wsmd(plot_data[[v]], plot_data$COUNTY_1, plot_data$iosw_lr))),
    stage = "Weighted (LR)"
  )
  
  smd_sl <- data.frame(
    variable = vars,
    smd = sapply(vars, function(v) abs(calc_wsmd(plot_data[[v]], plot_data$COUNTY_1, plot_data$iosw_sl))),
    stage = "Weighted (SL)"
  )
  
  smd_all <- bind_rows(smd_before, smd_lr, smd_sl) |>
    mutate(
      stage = factor(stage, levels = c("Unweighted", "Weighted (LR)", "Weighted (SL)")),
      variable = factor(variable, levels = rev(vars))
    )
  
  custom_colors <- c(
    "Unweighted"    = "#F8766D",
    "Weighted (LR)" = "#7570B3",
    "Weighted (SL)" = "#E69F00"
  )
  
  pd <- position_dodge(width = 0.4)
  
  # Generate Love plot
  p_smd <- ggplot(smd_all, aes(x = smd, y = variable, color = stage, shape = stage)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.7) +
    geom_point(size = 3.2, position = pd) +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = c("Unweighted" = 16, "Weighted (LR)" = 17, "Weighted (SL)" = 15)) +
    scale_y_discrete(labels = var_labels) +
    scale_x_continuous(
      limits = c(0, max(0.15, max(smd_all$smd, na.rm = TRUE) * 1.05)),
      breaks = seq(0, 1, by = 0.05)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray92", linetype = "dotted")
    ) +
    labs(
      x = "Absolute Standardized Mean Difference",
      y = NULL,
      color = "Weighting Method",
      shape = "Weighting Method"
    )
  
  return(p_smd)
}



#' Assess covariate overlap and transportability between source and target counties
#'
#' @description
#' Pipeline evaluating covariate balance when transporting estimates from source
#' counties to a target county. Computes propensity scores, Tipton's Beta 
#' overlap metrics, density graphs, and covariate balance Love plots.
#'
#' @param county_1 Character vector of source county names.
#' @param county_0 Character string or vector of target county name.
#' @param zsite.SL.library Character vector of SuperLearner algorithms for selection modeling. 
#'   Default is c('SL.glm', 'SL.ranger', 'SL.glmnet').
#' @param vars Character vector of environmental covariates. 
#'   Default is c("tmax", "prec", "hum", "shrtwv_rad", "wind").
#' @param var_labels Named character vector mapping raw variable names to new labels. 
#'  Default is c("tmax" = "Max. Temperature", "prec" = "Precipitation", 
#'  "hum" = "Relative Humidity", "shrtwv_rad" = "Shortwave Radiation",
#'  "wind" = "Wind Speed")

assess_transport_overlap <- function(county_1, 
                                     county_0, 
                                     zsite.SL.library = c('SL.glm', 'SL.ranger', 'SL.glmnet'),
                                     vars = c(
                                       "tmax", 
                                       "prec", 
                                       "hum", 
                                       "shrtwv_rad", 
                                       "wind"),
                                     var_labels = c(
                                       "tmax"       = "Max. Temperature",
                                       "prec"       = "Precipitation",
                                       "hum"        = "Relative Humidity",
                                       "shrtwv_rad" = "Shortwave Radiation",
                                       "wind"       = "Wind Speed"
                                     )) {
  
  # Prepare subset data
  county_1_data <- exposed_data |>
    mutate(COUNTY_1 = case_when(
      countyname %in% county_1 ~ 1,
      countyname %in% county_0 ~ 0,
      TRUE ~ NA_real_
    )) |>
    filter(!is.na(COUNTY_1))
  
  # Plot PS density overlap and calculate Tipton's Beta
  plot_data <- plot_ps_overlap(county_1_data, zsite.SL.library, vars)
  
  # Love plot
  plot_covariate_balance(plot_data, vars, var_labels)
}



#######################################################################################################
### Link between overlap and estimation error
#######################################################################################################

#' Evaluate the relationship between propensity score overlap and transport estimation error
#'
#' @description
#' Performs an iterative simulation across random combinations of target and source 
#' counties (from 1 up to max_sources) to evaluate how Tipton's Beta overlap
#' statistic relates to estimation errors (absolute and relative) of transported
#' effect estimates computed via Inverse Odds of Sampling Weights (IOSW).
#'
#' @param N Integer specifying the number of simulation iterations. Defaults to 200.
#' @param max_sources Integer specifying the maximum number of source counties to sample 
#'   jointly per iteration. Defaults to 5.
#' @param seed Integer specifying the random seed for reproducibility, or NULL. 
#'   Defaults to 42.
#' @param zsite.SL.library Character vector specifying SuperLearner algorithms for 
#'   selection propensity score modeling. Defaults to c('SL.glm', 'SL.ranger', 'SL.glmnet').
#'
#' @return A tibble with N rows containing simulation results:
#'   - iter : Iteration index ;
#'   - target : Character name of the target county ;
#'   - n_sources : Integer number of sampled source counties ;
#'   - beta : Tipton's Beta overlap statistic between combined sources and target ;
#'   - est_transported : Transported effect estimate (RR) via IOSW ;
#'   - est_observed : True observed effect estimate (RR) in the target county ;
#'   - abs_error : Absolute estimation error;
#'   - rel_error : Relative estimation error.

eval_overlap_error_link <- function(N = 200, 
                                    max_sources = 5, 
                                    seed = 42, 
                                    zsite.SL.library = c('SL.glm', 'SL.ranger', 'SL.glmnet')) {
  
  if (!is.null(seed)) set.seed(seed)
  
  counties_names <- unique(exposed_data$countyname)
  
  results_list <- vector("list", N)
  
  for (i in 1:N) {
    
    # Random selection of source (counties_1) and target (county_0) counties
    county_0 <- sample(counties_names, size = 1)
    n_sources <- sample(1:min(max_sources, length(counties_names) - 1), size = 1)
    counties_1 <- sample(setdiff(counties_names, county_0), size = n_sources)
    
    county_1_data <- exposed_data |>
      filter(countyname %in% c(county_0, counties_1)) |>
      mutate(COUNTY_1 = ifelse(countyname %in% counties_1, 1, 0))
    
    # Estimate selection propensity score using SuperLearner
    cps_sl <- SuperLearner(
      Y = county_1_data$COUNTY_1,
      X = county_1_data[, c("tmax", "prec", "hum", "shrtwv_rad", "wind")],
      SL.library = zsite.SL.library,
      family = 'binomial'
    )$SL.predict[, 1]
    
    p_county_1 <- mean(county_1_data$COUNTY_1 == 1)
    p_county_0 <- 1 - p_county_1
    
    county_1_data <- county_1_data |>
      mutate(iosw_sl = ifelse(COUNTY_1 == 1, (1 - cps_sl) / cps_sl * p_county_1 / p_county_0, 0))
    
    # Estimate transported effect (RR)
    outcome_mod_sl <- glm(
      resp ~ offset(log(total.pop)) + exp_day + time_unit + 
        tmax + prec + hum + shrtwv_rad + wind,
      family = quasipoisson,
      data = county_1_data,
      weights = iosw_sl,
      subset = (iosw_sl > 0)
    )
    
    estimate_sl <- exp(coef(outcome_mod_sl)["exp_day"])
    
    # Extract true observed effect for the target county
    observed <- exposed_county_effects |> 
      filter(county == county_0) |> 
      pull(estimate)
    
    beta <- calc_tipton_beta(cps_sl, county_1_data$COUNTY_1, common_support = FALSE)
    
    results_list[[i]] <- tibble(
      iter = i,
      target = county_0,
      n_sources = n_sources,
      beta = beta,
      est_transported = estimate_sl,
      est_observed = observed,
      abs_error = abs(estimate_sl - observed),
      rel_error = abs(estimate_sl - observed) / observed
    )
  }
  
  bind_rows(results_list)
}



#######################################################################################################
### Comparison between extrapolation and transportability
#######################################################################################################

#' Compute propensity score overlap between a target county and all candidate source counties
#'
#' @description
#' Calculates Tipton's Beta between a specified target county and all other available
#' candidate source counties using SuperLearner for selection propensity score estimation.
#'
#' @param county_0 Character string specifying the target county name.
#' @param zsite.SL.library Character vector specifying SuperLearner algorithms for 
#'   selection propensity score modeling. Defaults to c('SL.glm', 'SL.ranger', 'SL.glmnet').
#'
#' @return A data.frame with two columns:
#'   - county : Character name of the candidate source county ;
#'   - overlap : Tipton's Beta overlap statistic with county_0.

calc_source_overlaps <- function(county_0, zsite.SL.library = c('SL.glm', 'SL.ranger', 'SL.glmnet')) {
  
  counties_names <- setdiff(unique(exposed_data$countyname), county_0)
  
  overlap_results <- lapply(counties_names, function(county_1) {
    
    county_1_data <- exposed_data |>
      filter(countyname %in% c(county_0, county_1)) |>
      mutate(COUNTY_1 = if_else(countyname == county_1, 1, 0))
    
    # Estimate selection propensity score using SuperLearner
    cps_sl <- SuperLearner(
      Y = county_1_data$COUNTY_1,
      X = county_1_data[, c("tmax", "prec", "hum", "shrtwv_rad", "wind")],
      SL.library = zsite.SL.library,
      family = 'binomial'
    )$SL.predict[, 1]
    
    # Calculate Tipton's Beta
    overlap <- calc_tipton_beta(cps_sl, county_1_data$COUNTY_1, common_support = FALSE)
    
    data.frame(county = county_1, overlap = overlap)
  })
  
  bind_rows(overlap_results)
}


#' Evaluate transported and naive extrapolated effect estimates for a source selection scenario
#'
#' @description
#' Helper function that computes and formats point estimates and 95% confidence intervals 
#' for both transported effects and naive unweighted extrapolated effects for a given set
#' of source counties relative to a target county.
#'
#' @param source_counties Character vector of source county names included in the scenario.
#' @param label Character string specifying the descriptive scenario label.
#' @param county_0 Character string specifying the target county name.
#'
#' @return A data.frame containing combined estimation results with columns: Scenario, 
#' Method, estimate, lower_ci and upper_ci.

eval_transport_scenario <- function(source_counties, label, county_0) {
  
  # Transported effect estimate
  tr_raw <- transport_iosw(source_counties, county_0)[2, ]
  tr_df  <- data.frame(
    Scenario = label, Method = "Transportability",
    estimate = tr_raw$estimate, lower_ci = tr_raw$lower_ci, upper_ci = tr_raw$upper_ci
  )
  
  # Naive extrapolated effect estimate
  idx    <- which(exposed_data$countyname %in% source_counties)
  ex_raw <- estimate_effect(idx, count = TRUE)
  ex_df  <- data.frame(
    Scenario = label, Method = "Extrapolation",
    estimate = as.numeric(ex_raw["exp_day"]),
    lower_ci = as.numeric(ex_raw["2.5 %"]),
    upper_ci = as.numeric(ex_raw["97.5 %"])
  )
  
  bind_rows(tr_df, ex_df)
}


#' Compare naive extrapolation and causal transportability across source selection scenarios
#'
#' @description
#' Evaluates and visualizes the performance of naive extrapolation vs causal 
#' transportability across 6 distinct source selection strategies (all counties, 
#' 5 most/least similar by Tipton's Beta overlap, 5 random, and 5 highest/lowest 
#' observed effect estimates) against the target county's true observed effect.
#'
#' @param county_0 Character string specifying the target county name.
#' @param seed Integer specifying the random seed for sampling reproducibility, or NULL. 
#'   Defaults to 123.
#'
#' @return A ggplot object displaying side-by-side point estimates and 95\% 
#'  confidence intervals across scenarios and estimation methods.

plot_extrapolation_vs_transport <- function(county_0, seed = 123) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Extract true observed effect in target county
  true_val <- filter(exposed_county_effects, county == county_0)
  
  # Calculate propensity score overlap for all candidate source counties
  overlaps_df <- calc_source_overlaps(county_0)
  
  # Define candidate source populations
  all_counties   <- unique(exposed_data$countyname)
  other_counties <- setdiff(all_counties, county_0)
  
  n_select <- min(5, nrow(overlaps_df))
  least_similar_5 <- overlaps_df |> 
    arrange(overlap) |> 
    slice_head(n = n_select) |> 
    pull(county)
  most_similar_5 <- overlaps_df |> 
    arrange(desc(overlap)) |> 
    slice_head(n = n_select) |> 
    pull(county)
  
  random_5 <- sample(other_counties, size = min(5, length(other_counties)))
  
  # Construct source selection scenarios
  scenarios <- list(
    "All other counties"     = other_counties,
    "5 most similar"         = most_similar_5,
    "5 least similar"        = least_similar_5,
    "5 random counties"      = random_5,
    "5 highest effect size"  = top_5_counties,
    "5 lowest effect size"   = no_effect_counties
  )
  
  # Evaluate all scenarios
  results <- bind_rows(lapply(names(scenarios), function(lbl) {
    eval_transport_scenario(scenarios[[lbl]], lbl, county_0 = county_0)
  }))
  
  results$Scenario <- factor(
    results$Scenario,
    levels = c("5 least similar", 
               "5 random counties", 
               "All other counties", 
               "5 most similar",
               "5 highest effect size",
               "5 lowest effect size")
  )
  
  # Generate comparison plot
  pd <- position_dodge(width = 0.5)
  
  ggplot(results, aes(x = estimate, y = Scenario, color = Method)) +
    annotate("rect", 
             xmin = true_val$lower_ci, xmax = true_val$upper_ci, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "gray20") +
    geom_vline(xintercept = true_val$estimate, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_errorbar(
      aes(xmin = lower_ci, xmax = upper_ci), 
      position = pd, 
      height = 0.25, 
      linewidth = 0.8
    ) +
    geom_point(
      aes(shape = Method), 
      position = pd, 
      size = 2.5
    ) +
    scale_color_manual(values = c("Extrapolation" = "#d73027", "Transportability" = "#4575b4")) +
    scale_x_continuous(breaks = seq(-3, 3, by = 0.1)) +
    labs(
      title = paste("Estimated Causal Effect (RR) — Target:", county_0),
      subtitle = "Dashed line (shaded area) = Target county observed effect (95% CI)",
      x = "Estimated Rate Ratio (RR)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major.y = element_blank()
    )
}