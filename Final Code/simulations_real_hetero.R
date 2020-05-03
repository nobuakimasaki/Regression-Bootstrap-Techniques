library(tidyverse)
library(parallel)
library(boot)

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

# sampling cases - BCA
cases_bca <- function(y_sample, x_sample, model, n, m){
  stats <- function(xy_sample, i) {
    xstar <- xy_sample[i, 1]
    ystar <- xy_sample[i, 2]
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
  }
  boot.out <- boot(cbind(x_sample, y_sample), statistic=stats, R = m)
  
  itct_lower_cases <- boot.ci(boot.out = boot.out, type = "bca", index = 1)$bca[c(4)]
  itct_upper_cases <- boot.ci(boot.out = boot.out, type = "bca", index = 1)$bca[c(5)]
  slope_lower_cases <- boot.ci(boot.out = boot.out, type = "bca", index = 2)$bca[c(4)]
  slope_upper_cases <- boot.ci(boot.out = boot.out, type = "bca", index = 2)$bca[c(5)]
  sigma_lower_cases <- boot.ci(boot.out = boot.out, type = "bca", index = 3)$bca[c(4)]
  sigma_upper_cases <- boot.ci(boot.out = boot.out, type = "bca", index = 3)$bca[c(5)]
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

# sampling residuals - BCA
boot_resid_bca <- function(y_sample, x_sample, model, n, m){
  regstats <- function(dat, i) {
    #dat is a data frame (r, x, yhat)
    #r are the modified centered residuals, yhat are the fits
    ystar <- dat$yhat + dat$r[i]
    xstar <- dat$x
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope = Lb$coeff[2], s = s)
  }
  
  m.resid <- rstandard(model, sd = 1)
  r <- m.resid - mean(m.resid)
  yhat <- model$fitted.values
  dat <- data.frame(r = r, x = x_sample, yhat = yhat)
  boot.out <- boot(dat, statistic=regstats, R = m)
  
  itct_lower_boot_resid <- boot.ci(boot.out = boot.out, type = "bca", index = 1)$bca[c(4)]
  itct_upper_boot_resid <- boot.ci(boot.out = boot.out, type = "bca", index = 1)$bca[c(5)]
  slope_lower_boot_resid <- boot.ci(boot.out = boot.out, type = "bca", index = 2)$bca[c(4)]
  slope_upper_boot_resid <- boot.ci(boot.out = boot.out, type = "bca", index = 2)$bca[c(5)]
  sigma_lower_boot_resid <- boot.ci(boot.out = boot.out, type = "bca", index = 3)$bca[c(4)]
  sigma_upper_boot_resid <- boot.ci(boot.out = boot.out, type = "bca", index = 3)$bca[c(5)]
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

# wild - BCA
wild_bca <- function(y_sample, x_sample, model, n, m){
  regstats <- function(dat, i) {
    #dat is a data frame (r, x, yhat)
    #r are the residuals, yhat are the fits
    e <- dat$r
    token <- rbinom(n = n, size = 1, p = (5 + sqrt(5))/10)
    estar <- ifelse(token == 1, e*(1-sqrt(5)/2), e*(1+sqrt(5)/2))
    ystar <- dat$yhat + estar
    xstar <- dat$x
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope = Lb$coeff[2], s = s)
  }
  
  r <- model$residuals
  yhat <- model$fitted.values
  dat <- data.frame(r = r, x = x_sample, yhat = yhat)
  boot.out <- boot(dat, statistic=regstats, R = m)  # m should be big
  
  itct_lower_wild <- boot.ci(boot.out = boot.out, type = "bca", index = 1)$bca[c(4)]
  itct_upper_wild <- boot.ci(boot.out = boot.out, type = "bca", index = 1)$bca[c(5)]
  slope_lower_wild <- boot.ci(boot.out = boot.out, type = "bca", index = 2)$bca[c(4)]
  slope_upper_wild <- boot.ci(boot.out = boot.out, type = "bca", index = 2)$bca[c(5)]
  sigma_lower_wild <- boot.ci(boot.out = boot.out, type = "bca", index = 3)$bca[c(4)]
  sigma_upper_wild <- boot.ci(boot.out = boot.out, type = "bca", index = 3)$bca[c(5)]
  return(c(itct_lower_wild, itct_upper_wild, slope_lower_wild, slope_upper_wild,
           sigma_lower_wild, sigma_upper_wild))
}


runFunc <- function(f, y_sample, x_sample, model, n, m){
  if (f == 1) {return(cases_conf_int(y_sample, x_sample, model, n, m))}
  if (f == 2) {return(boot_resid_conf_int(y_sample, x_sample, model, n, m))}
  if (f == 3) {return(wild_conf_int(y_sample, x_sample, model, n, m))}
  if (f == 4) {return(cases_bca(y_sample, x_sample, model, n, m))}
  if (f == 5) {return(boot_resid_bca(y_sample, x_sample, model, n, m))}
  if (f == 6) {return(wild_bca(y_sample, x_sample, model, n, m))}
}

sim_iter <- function(d_, y, x) {
  
  slope <- d_$slope
  intercept <- d_$intercept
  sigma <- d_$sigma
  n <- d_$n
  # sampling from population
  i <- sample(1:N, replace = FALSE, size = n)
  y_sample <- y[i]
  x_sample <- x[i]
  model <- lm(y_sample ~ x_sample)
  
  conf_limits_model <- classical_conf_int(y_sample, x_sample, model)
  itct_lower_model <- conf_limits_model[1]
  itct_upper_model <- conf_limits_model[2]
  slope_lower_model <- conf_limits_model[3]
  slope_upper_model <- conf_limits_model[4]
  
  if (itct_lower_model > intercept) { # interval is too high
    result_itct_model <- 1
  } else if (intercept > itct_upper_model) {
    result_itct_model <- -1
  } else {result_itct_model <- 0}
  
  if (slope_lower_model > slope) {
    result_slope_model <- 1
  } else if (slope > slope_upper_model) {
    result_slope_model <- -1
  } else {result_slope_model <- 0}
  
  results <- numeric(0)
  for (f in 1:6){
    confs <- runFunc(1, y_sample, x_sample, model, n, m)
    itct_lower<- confs[1]
    itct_upper <- confs[2]
    slope_lower <- confs[3]
    slope_upper <- confs[4]
    sigma_lower <- confs[5]
    sigma_upper <- confs[6]
    
    if (itct_lower > intercept) {  # interval too high
      result_itct <- 1
    } else if (intercept > itct_upper) {
      result_itct <- -1
    } else {result_itct <- 0}
    
    if (slope_lower > slope) {
      result_slope <- 1
    } else if (slope > slope_upper) {
      result_slope <- -1
    } else {result_slope <- 0}
    
    if (sigma_lower > sigma) {
      result_sigma <- 1
    } else if (sigma > sigma_upper) {
      result_sigma <- -1
    } else {result_sigma <- 0}
    
    results <- c(results, result_itct, result_slope, result_sigma)
  }
  
  print("finished 1 iteration")
  
  return(c(result_itct_model, result_slope_model, results))
}

get_coverage <- function(x) {
  coverage <- mean(x == 0)
  too_high <- mean(x == 1)
  too_low <- mean(x == -1)
  return(c(coverage, too_high, too_low))
}

set.seed(17)

# loading data
ASIE_data <- read.csv("ASIE2013.csv")

# creating population
ASIE_data$log_total_assets <- log(ASIE_data$total_assets)
ASIE_data$net_income_scale <- scale(ASIE_data$net_income)

ASIE_data <- ASIE_data %>% filter(net_income_scale < 38 & log_total_assets > 12.5)

ggplot(ASIE_data, aes(x = log_total_assets, y = net_income_scale)) + 
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Log Total Assets", y = "Net Income (Normalized)") 

x <- ASIE_data$log_total_assets
y <- ASIE_data$net_income_scale

iter <- 500
N <- length(x)
sample_size <- c(60, 120, 240)
m <- 1000

pop_model <- lm(y ~ x)
slope <- coef(pop_model)[2]
intercept <- coef(pop_model)[1]
sigma <- summary(pop_model)$sigma

for (n in sample_size) {
  
  df <- data.frame("slope" = slope, "intercept" = intercept, "sigma" = sigma, "n" = n)
  df <- df[rep(seq_len(nrow(df)), each = iter), ]
  
  d_ <- list()
  for (i in seq(nrow(df))) {
    l <- as.list(df[i,])
    d_[[i]] <- l
  }
  
  sim <- mclapply(d_, sim_iter, y = y, x = x, mc.cores = 1)
  sim_df <- as.data.frame(t(matrix(unlist(sim), nrow=length(unlist(sim[1])))))
  sim_df_summary <- as.data.frame(apply(sim_df, 2, get_coverage))
  colnames(sim_df_summary) <- c("itct_model", "slope_model", 
                                "itct_cases", "slope_cases", "sigma_cases",
                                "itct_cases_bca", "slope_cases_bca", "sigma_cases_bca",
                                "itct_boot_resid", "slope_boot_resid", "sigma_boot_resid",
                                "itct_boot_resid_bca", "slope_boot_resid_bca", "sigma_boot_resid_bca",
                                "itct_wild", "slope_wild", "sigma_wild",
                                "itct_wild_bca", "slope_wild_bca", "sigma_wild_bca")
  rownames(sim_df_summary) <- c("coverage", "too_high", "too_low")
  
  file_name <- paste("result_real_hetero_data_", n, ".csv", sep = "")
  
  write.csv(sim_df_summary, file_name)
}
