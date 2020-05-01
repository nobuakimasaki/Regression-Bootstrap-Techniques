library(tidyverse)
library(parallel)

# standard confidence interval
classical_conf_int <- function(y_sample, x_sample, model){
  
  itct_lower_model <- confint(model)[[1]]
  itct_upper_model <- confint(model)[[3]]
  slope_lower_model <- confint(model)[[2]]
  slope_upper_model <- confint(model)[[4]]
  return(c(itct_lower_model, itct_upper_model, slope_lower_model, slope_upper_model))
}

# sampling cases
cases_conf_int <- function(y_sample, x_sample, model, n, m){
  
  out_cases <- replicate(m, expr = {
    j <- sample(1:n, replace = TRUE, size = n)
    xstar <- x_sample[j]
    ystar <- y_sample[j]
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope = Lb$coeff[2], s = s)
  })
  
  itct_lower_cases <- quantile(out_cases[1,], probs = 0.025)
  itct_upper_cases <- quantile(out_cases[1,], probs = 0.975)
  slope_lower_cases <- quantile(out_cases[2,], probs = 0.025)
  slope_upper_cases <- quantile(out_cases[2,], probs = 0.975)
  sigma_lower_cases <- quantile(out_cases[3,], probs = 0.025)
  sigma_upper_cases <- quantile(out_cases[3,], probs = 0.975)
  return(c(itct_lower_cases, itct_upper_cases, slope_lower_cases, slope_upper_cases,
           sigma_lower_cases, sigma_upper_cases))
}

# sampling residuals
boot_resid_conf_int <- function(y_sample, x_sample, model, n, m){
  
  m.resid <- rstandard(model, sd = 1)
  r <- m.resid - mean(m.resid)
  out_boot_resid <- replicate(m, expr = {
    estar <- sample(r, replace = TRUE, size = n)
    ystar <- model$fitted.values + estar
    Lb <- lm(ystar ~ x_sample)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope = Lb$coeff[2], s = s)
  })
  
  itct_lower_boot_resid <- quantile(out_boot_resid[1,], probs = 0.025)
  itct_upper_boot_resid <- quantile(out_boot_resid[1,], probs = 0.975)
  slope_lower_boot_resid <- quantile(out_boot_resid[2,], probs = 0.025)
  slope_upper_boot_resid <- quantile(out_boot_resid[2,], probs = 0.975)
  sigma_lower_boot_resid <- quantile(out_boot_resid[3,], probs = 0.025)
  sigma_upper_boot_resid <- quantile(out_boot_resid[3,], probs = 0.975)
  return(c(itct_lower_boot_resid, itct_upper_boot_resid, slope_lower_boot_resid, slope_upper_boot_resid,
           sigma_lower_boot_resid, sigma_upper_boot_resid))
}

# wild
wild_conf_int <- function(y_sample, x_sample, model, n, m){
  
  out_wild <- replicate(m, expr = {
    e <- model$residuals
    token <- rbinom(n = n, size = 1, p = (5 + sqrt(5))/10)
    estar <- ifelse(token == 1, e*(1-sqrt(5)/2), e*(1+sqrt(5)/2))
    ystar <- model$fitted.values + estar
    Lb <- lm(ystar ~ x_sample)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope = Lb$coeff[2], s = s)
  })
  itct_lower_wild <- quantile(out_wild[1,], probs = 0.025)
  itct_upper_wild <- quantile(out_wild[1,], probs = 0.975)
  slope_lower_wild <- quantile(out_wild[2,], probs = 0.025)
  slope_upper_wild <- quantile(out_wild[2,], probs = 0.975)
  sigma_lower_wild <- quantile(out_wild[3,], probs = 0.025)
  sigma_upper_wild <- quantile(out_wild[3,], probs = 0.975)
  return(c(itct_lower_wild, itct_upper_wild, slope_lower_wild, slope_upper_wild,
           sigma_lower_wild, sigma_upper_wild))
}

sim_iter <- function(d_, y, x, outlier) {
  
  slope <- d_$slope
  intercept <- d_$intercept
  sigma <- d_$sigma
  n <- d_$n
  
  # sampling from population
  y_without_outlier <- y[-outlier]
  x_without_outlier <- x[-outlier]
  indices <- 1:N
  avail_indices <- indices[-outlier]
  i <- sample(avail_indices, replace = FALSE, size = 232)
  y_sample <- c(y_without_outlier[i], y[outlier])
  x_sample <- c(x_without_outlier[i], x[outlier])
  model <- lm(y_sample ~ x_sample)
  
  conf_limits_model <- classical_conf_int(y_sample, x_sample, model)
  itct_lower_model <- conf_limits_model[1]
  itct_upper_model <- conf_limits_model[2]
  slope_lower_model <- conf_limits_model[3]
  slope_upper_model <- conf_limits_model[4]
  
  if (itct_lower_model > intercept) {result_itct_model <- 1} # interval is too high
  else if (intercept > itct_upper_model) {result_itct_model <- -1}
  else {result_itct_model <- 0}
  if (slope_lower_model > slope) {result_slope_model <- 1}
  else if (slope > slope_upper_model) {result_slope_model <- -1}
  else {result_slope_model <- 0}
  
  conf_limits_cases <- cases_conf_int(y_sample, x_sample, model, n, m)
  itct_lower_cases <- conf_limits_cases[1]
  itct_upper_cases <- conf_limits_cases[2]
  slope_lower_cases <- conf_limits_cases[3]
  slope_upper_cases <- conf_limits_cases[4]
  sigma_lower_cases <- conf_limits_cases[5]
  sigma_upper_cases <- conf_limits_cases[6]
  
  if (itct_lower_cases > intercept) {result_itct_cases <- 1} # interval is too high
  else if (intercept > itct_upper_cases) {result_itct_cases <- -1}
  else {result_itct_cases <- 0}
  if (slope_lower_cases > slope) {result_slope_cases <- 1}
  else if (slope > slope_upper_cases) {result_slope_cases <- -1}
  else {result_slope_cases <- 0}
  if (sigma_lower_cases > sigma) {result_sigma_cases <- 1}
  else if (sigma > sigma_upper_cases) {result_sigma_cases <- -1}
  else {result_sigma_cases <- 0}
  
  conf_limits_boot_resid <- boot_resid_conf_int(y_sample, x_sample, model, n, m)
  itct_lower_boot_resid <- conf_limits_boot_resid[1]
  itct_upper_boot_resid <- conf_limits_boot_resid[2]
  slope_lower_boot_resid <- conf_limits_boot_resid[3]
  slope_upper_boot_resid <- conf_limits_boot_resid[4]
  sigma_lower_boot_resid <- conf_limits_boot_resid[5]
  sigma_upper_boot_resid <- conf_limits_boot_resid[6]
  
  if (itct_lower_boot_resid > intercept) {result_itct_boot_resid <- 1} # interval is too high
  else if (intercept > itct_upper_boot_resid) {result_itct_boot_resid <- -1}
  else {result_itct_boot_resid <- 0}
  if (slope_lower_boot_resid > slope) {result_slope_boot_resid <- 1}
  else if (slope > slope_upper_boot_resid) {result_slope_boot_resid <- -1}
  else {result_slope_boot_resid <- 0}
  if (sigma_lower_boot_resid > sigma) {result_sigma_boot_resid <- 1}
  else if (sigma > sigma_upper_boot_resid) {result_sigma_boot_resid <- -1}
  else {result_sigma_boot_resid <- 0}
  
  conf_limits_wild <- wild_conf_int(y_sample, x_sample, model, n, m)
  itct_lower_wild <- conf_limits_wild[1]
  itct_upper_wild <- conf_limits_wild[2]
  slope_lower_wild <- conf_limits_wild[3]
  slope_upper_wild <- conf_limits_wild[4]
  sigma_lower_wild <- conf_limits_wild[5]
  sigma_upper_wild <- conf_limits_wild[6]
  
  if (itct_lower_wild > intercept) {result_itct_wild <- 1} # interval is too high
  else if (intercept > itct_upper_wild) {result_itct_wild <- -1}
  else {result_itct_wild <- 0}
  if (slope_lower_wild > slope) {result_slope_wild <- 1}
  else if (slope > slope_upper_wild) {result_slope_wild <- -1}
  else {result_slope_wild <- 0}
  if (sigma_lower_wild > sigma) {result_sigma_wild <- 1}
  else if (sigma > sigma_upper_wild) {result_sigma_wild <- -1}
  else {result_sigma_wild <- 0}
  
  print("finished 1 iteration")
  
  return(c(result_itct_model, result_slope_model, 
           result_itct_cases, result_slope_cases, result_sigma_cases,
           result_itct_boot_resid, result_slope_boot_resid, result_sigma_boot_resid,
           result_itct_wild, result_slope_wild, result_sigma_wild))
}

get_coverage <- function(x) {
  coverage <- mean(x == 0)
  too_high <- mean(x == 1)
  too_low <- mean(x == -1)
  return(c(coverage, too_high, too_low))
}

set.seed(17)

# loading data
NIDS_data <- read_csv("NIDS_data.csv") %>%
  filter(w5_best_age_yrs > 0) %>%
  filter(log(w5_pi_tot_ass_i) > 0)

x <- NIDS_data$w5_best_age_yrs
y <- log(NIDS_data$w5_pi_tot_ass_i)

sort(y)[(length(y)-3):length(y)]
which(y > 17.90290) # 256 1891 2615 3686
y[c(256,3686,2615,1891)]
x[c(256,3686,2615,1891)]

sort(y)[1:4]
which(y < 1.386295) # 2969  4609  7803 10040
x[c(4609,2969,7803,10040)]
y[c(4609,2969,7803,10040)]

iter <- 500
N <- length(x)
sample_size <- 240
m <- 10000

# fitting lm and obtaining parameter values
pop_model <- lm(y ~ x)
slope <-  coef(pop_model)[2]
intercept <- coef(pop_model)[1]
sigma <- summary(pop_model)$sigma

# example of a sample
outlier <- c(256,3686,2615,1891,4609,2969,7803,10040)
y_without_outlier <- y[-outlier]
x_without_outlier <- x[-outlier]
indices <- 1:N
avail_indices <- indices[-outlier]
i <- sample(avail_indices, replace = FALSE, size = sample_size[1]-8)
y_sample <- c(y_without_outlier[i], y[outlier])
x_sample <- c(x_without_outlier[i], x[outlier])
plot(y_sample ~ x_sample, col="lightblue", pch=19, cex=2,)
text(y_sample ~ x_sample, labels=1:sample_size[1], cex=0.8, font=1)

for (n in sample_size) {
  
  df <- data.frame("slope" = slope, "intercept" = intercept, "sigma" = sigma, "n" = n)
  df <- df[rep(seq_len(nrow(df)), each = iter), ]
  
  d_ <- list()
  for (i in seq(nrow(df))) {
    l <- as.list(df[i,])
    d_[[i]] <- l
  }
  
  sim <- mclapply(d_, sim_iter, y = y, x = x, outlier = outlier, mc.cores = 1)
  sim_df <- as.data.frame(t(matrix(unlist(sim), nrow=length(unlist(sim[1])))))
  sim_df_summary <- as.data.frame(apply(sim_df, 2, get_coverage))
  colnames(sim_df_summary) <- c("itct_model", "slope_model", 
                                "itct_cases", "slope_cases", "sigma_cases",
                                "itct_boot_resid", "slope_boot_resid", "sigma_boot_resid",
                                "itct_wild", "slope_wild", "sigma_wild")
  rownames(sim_df_summary) <- c("coverage", "too_high", "too_low")
  
  file_name <- paste("result_outlier_", n, ".csv", sep = "")
  
  write.csv(sim_df_summary, file_name)
}
