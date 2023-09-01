library(readr)
library(gTests)
library(depthtestr)
library(vegan)
library(mvtnorm)

N1 = 50
N2 = 50
nsamp <- 500

alpha = 0.50
mu = 2
d = 10

thetas <- seq(0.25, 1, 0.125)

type_1_error_hd <- c()
type_1_error_mah <- c()
type_1_error_lp <- c()
type_1_error_pd <- c()
type_1_error_lcd_variant <- c()
type_1_error_lcd <- c()
type_1_error_cf <- c()
type_1_error_pvb <- c()

power_hd <- c()
power_mah <- c()
power_lp <- c()
power_pd <- c()
power_lcd_variant <- c()
power_lcd <- c()
power_cf <- c()
power_pvb <- c()

power_adjusted_hd <- c()
power_adjusted_mah <- c()
power_adjusted_pvb <- c()
power_adjusted_lp <- c()
power_adjusted_lcd_variant <- c()
power_adjusted_lcd <- c()
power_adjusted_pd <- c()

# type I error control
delta <- 0
for(j in 1:length(thetas)){
  
  set.seed(j)
  theta <- thetas[j]
  
  pvals_hd <- c()
  pvals_mah <- c()
  pvals_lp <- c()
  pvals_pd <- c()
  pvals_lcd_variant <- c()
  pvals_lcd <- c()
  pvals_pvb <- c()
  pvals_cf <- c()
  
  for(i in 1:nsamp){
    
    samp1_binom <- rbinom(N1, 1, alpha)
    samp2_binom <- rbinom(N2, 1, alpha)
    
    samp1 = rmvnorm(N1, mean = rep(0, d))*matrix(rep(samp1_binom, d), nrow=N1) +
      rmvnorm(N1, mean = rep(mu, d),
              sigma = diag(rep(theta^2, d)))*matrix(rep((1-samp1_binom), d), nrow=N1)
    samp2 = rmvnorm(N2, mean = rep(0, d),
                    sigma=diag(rep((1 + delta/sqrt(N1 + N2))^2, d)))*matrix(rep(samp2_binom, d), nrow=N2) +
      rmvnorm(N2, mean = rep(mu, d),
              sigma=diag(rep((theta + delta/sqrt(N1 + N2))^2, d)))*matrix(rep((1-samp2_binom), d), nrow=N2)
    
    
    pvals_hd[i] <- suppressMessages(depthTest(samp1, samp2, "halfspace",
                                              loo_correction = F)$pval)
    pvals_mah[i] <- suppressMessages(depthTest(samp1, samp2, "mahalanobis",
                                               loo_correction = F)$pval)
    pvals_lp[i] <- suppressMessages(depthTest(samp1, samp2, "lp",
                                              loo_correction = F)$pval)
    pvals_pd[i] <- suppressMessages(depthTest(samp1, samp2, "pd",
                                              loo_correction = F)$pval)
    pvals_lcd_variant[i] <- suppressMessages(depthTest(samp1, samp2, "lcd-variant",
                                               loo_correction = F)$pval)
    pvals_lcd[i] <- suppressMessages(depthTest(samp1, samp2, "lcd",
                                                loo_correction = F)$pval)
    pvals_pvb[i] <- suppressMessages(depthTest(samp1, samp2, "pvb",
                                               loo_correction = F)$pval)
    
    points <- rbind(samp1, samp2)
    x <- as.matrix(vegdist(points, binary=FALSE, method="euclidean"))
    child <- spantree(x)$kid
    edges <- t(rbind(c(2:(N1+N2)), child))
    
    pvals_cf[i] <- g.tests(edges, c(1:N1), c((N1+1):(N1+N2)), test.type="g", perm=0)$generalized$pval.approx
    
    
    print(i)
  }
  
  
  type_1_error_hd[j] <- mean(pvals_hd < 0.05)
  type_1_error_mah[j] <- mean(pvals_mah < 0.05)
  type_1_error_lp[j] <- mean(pvals_lp < 0.05)
  type_1_error_pd[j] <- mean(pvals_pd < 0.05)
  type_1_error_lcd_variant[j] <- mean(pvals_lcd_variant < 0.05)
  type_1_error_lcd[j] <- mean(pvals_lcd < 0.05)
  type_1_error_cf[j] <- mean(pvals_cf < 0.05)
  type_1_error_pvb[j] <- mean(pvals_pvb < 0.05)
  
  print(theta)
}

# power
delta <- 2
for(j in 1:length(thetas)){
  
  set.seed(j)
  theta <- thetas[j]
  
  pvals_hd <- c()
  pvals_mah <- c()
  pvals_lp <- c()
  pvals_pd <- c()
  pvals_lcd_variant <- c()
  pvals_lcd <- c()
  pvals_pvb <- c()
  pvals_cf <- c()
  
  stats_hd <- c()
  stats_null_hd <- c()
  stats_mah <- c()
  stats_null_mah <- c()
  stats_pvb <- c()
  stats_null_pvb <- c()
  stats_lp <- c()
  stats_null_lp <- c()
  stats_lcd_variant <- c()
  stats_null_lcd_variant <- c()
  stats_lcd <- c()
  stats_null_lcd <- c()
  stats_pd <- c()
  stats_null_pd <- c()
  
  for(i in 1:nsamp){
    
    samp1_binom <- rbinom(N1, 1, alpha)
    samp2_binom <- rbinom(N2, 1, alpha)
    samp3_binom <- rbinom(N2, 1, alpha)
    
    samp1 = rmvnorm(N1, mean = rep(0, d))*matrix(rep(samp1_binom, d), nrow=N1) +
      rmvnorm(N1, mean = rep(mu, d),
              sigma = diag(rep(theta^2, d)))*matrix(rep((1-samp1_binom), d), nrow=N1)
    samp2 = rmvnorm(N2, mean = rep(0, d),
                    sigma=diag(rep((1 + delta/sqrt(N1 + N2))^2, d)))*matrix(rep(samp2_binom, d), nrow=N2) +
      rmvnorm(N2, mean = rep(mu, d),
              sigma=diag(rep((theta + delta/sqrt(N1 + N2))^2, d)))*matrix(rep((1-samp2_binom), d), nrow=N2)
    samp3 = rmvnorm(N2, mean = rep(0, d))*matrix(rep(samp3_binom, d), nrow=N2) +
      rmvnorm(N2, mean = rep(mu, d),
              sigma = diag(rep(theta^2, d)))*matrix(rep((1-samp3_binom), d), nrow=N2)
    
    test_hd_12 <- suppressMessages(depthTest(samp1, samp2, "halfspace",
                                             loo_correction = F))
    test_mah_12 <- suppressMessages(depthTest(samp1, samp2, "mahalanobis",
                                              loo_correction = F))
    test_pvb_12 <- suppressMessages(depthTest(samp1, samp2, "pvb",
                                              loo_correction = F))
    test_lp_12 <- suppressMessages(depthTest(samp1, samp2, "lp",
                                             loo_correction = F))
    test_lcd_variant_12 <- suppressMessages(depthTest(samp1, samp2, "lcd-variant",
                                              loo_correction = F))
    test_lcd_12 <- suppressMessages(depthTest(samp1, samp2, "lcd",
                                               loo_correction = F))
    test_pd_12 <- suppressMessages(depthTest(samp1, samp2, "pd",
                                             loo_correction = F))
    
    test_hd_13 <- suppressMessages(depthTest(samp1, samp3, "halfspace",
                                             loo_correction = F))
    test_mah_13 <- suppressMessages(depthTest(samp1, samp3, "mahalanobis",
                                              loo_correction = F))
    test_pvb_13 <- suppressMessages(depthTest(samp1, samp3, "pvb",
                                              loo_correction = F))
    test_lp_13 <- suppressMessages(depthTest(samp1, samp3, "lp",
                                             loo_correction = F))
    test_lcd_variant_13 <- suppressMessages(depthTest(samp1, samp3, "lcd-variant",
                                              loo_correction = F))
    test_lcd_13 <- suppressMessages(depthTest(samp1, samp3, "lcd",
                                               loo_correction = F))
    test_pd_13 <- suppressMessages(depthTest(samp1, samp3, "pd",
                                             loo_correction = F))
    
    stats_hd[i] <- test_hd_12$teststat
    stats_null_hd[i] <- test_hd_13$teststat
    stats_mah[i] <- test_mah_12$teststat
    stats_null_mah[i] <- test_mah_13$teststat
    stats_pvb[i] <- test_pvb_12$teststat
    stats_null_pvb[i] <- test_pvb_13$teststat
    stats_lp[i] <- test_lp_12$teststat
    stats_null_lp[i] <- test_lp_13$teststat
    stats_lcd_variant[i] <- test_lcd_variant_12$teststat
    stats_null_lcd_variant[i] <- test_lcd_variant_13$teststat
    stats_lcd[i] <- test_lcd_12$teststat
    stats_null_lcd[i] <- test_lcd_13$teststat
    stats_pd[i] <- test_pd_12$teststat
    stats_null_pd[i] <- test_pd_13$teststat
    
    pvals_hd[i] <- test_hd_12$pval
    pvals_mah[i] <- test_mah_12$pval
    pvals_lp[i] <- test_lp_12$pval
    pvals_pd[i] <- test_pd_12$pval
    pvals_lcd_variant[i] <- test_lcd_variant_12$pval
    pvals_lcd[i] <- test_lcd_12$pval
    pvals_pvb[i] <- test_pvb_12$pval
    
    points <- rbind(samp1, samp2)
    x <- as.matrix(vegdist(points, binary=FALSE, method="euclidean"))
    child <- spantree(x)$kid
    edges <- t(rbind(c(2:(N1+N2)), child))
    
    pvals_cf[i] <- g.tests(edges, c(1:N1), c((N1+1):(N1+N2)), test.type="g", perm=0)$generalized$pval.approx
    
    
    
    print(i)
  }
  
  
  power_hd[j] <- mean(pvals_hd < 0.05)
  power_mah[j] <- mean(pvals_mah < 0.05)
  power_lp[j] <- mean(pvals_lp < 0.05)
  power_pd[j] <- mean(pvals_pd < 0.05)
  power_lcd_variant[j] <- mean(pvals_lcd_variant < 0.05)
  power_lcd[j] <- mean(pvals_lcd < 0.05)
  power_cf[j] <- mean(pvals_cf < 0.05)
  power_pvb[j] <- mean(pvals_pvb < 0.05)
  
  power_adjusted_hd[j] <- mean(stats_hd > quantile(stats_null_hd, 0.975)) +
    mean(stats_hd < quantile(stats_null_hd, 0.025))
  power_adjusted_mah[j] <- mean(stats_mah > quantile(stats_null_mah, 0.975)) +
    mean(stats_mah < quantile(stats_null_mah, 0.025))
  power_adjusted_pvb[j] <- mean(stats_pvb > quantile(stats_null_pvb, 0.975)) +
    mean(stats_pvb < quantile(stats_null_pvb, 0.025))
  power_adjusted_lp[j] <- mean(stats_lp > quantile(stats_null_lp, 0.975)) +
    mean(stats_lp < quantile(stats_null_lp, 0.025))
  power_adjusted_pd[j] <- mean(stats_pd > quantile(stats_null_pd, 0.975)) +
    mean(stats_pd < quantile(stats_null_pd, 0.025))
  power_adjusted_lcd_variant[j] <- mean(stats_lcd_variant > quantile(stats_null_lcd_variant, 0.975)) +
    mean(stats_lcd_variant < quantile(stats_null_lcd_variant, 0.025))
  power_adjusted_lcd[j] <- mean(stats_lcd > quantile(stats_null_lcd, 0.975)) +
    mean(stats_lcd < quantile(stats_null_lcd, 0.025))
  
  print(theta)
}


output <- data.frame(theta = thetas,
                     type_1_error_hd = type_1_error_hd,
                     type_1_error_mah = type_1_error_mah,
                     type_1_error_lp = type_1_error_lp,
                     type_1_error_pd = type_1_error_pd,
                     type_1_error_lcd_variant = type_1_error_lcd_variant,
                     type_1_error_lcd = type_1_error_lcd,
                     type_1_error_pvb = type_1_error_pvb,
                     type_1_error_cf = type_1_error_cf,
                     power_hd = power_hd,
                     power_mah = power_mah,
                     power_lp = power_lp,
                     power_pd = power_pd,
                     power_lcd_variant = power_lcd_variant,
                     power_lcd = power_lcd,
                     power_pvb = power_pvb,
                     power_cf = power_cf,
                     power_adjusted_hd = power_adjusted_hd,
                     power_adjusted_mah = power_adjusted_mah,
                     power_adjusted_lp = power_adjusted_lp,
                     power_adjusted_pd = power_adjusted_pd,
                     power_adjusted_lcd_variant = power_adjusted_lcd_variant,
                     power_adjusted_lcd = power_adjusted_lcd,
                     power_adjusted_pvb = power_adjusted_pvb)

write_csv(output, file = "mixture_results_power_vs_sd_uncorrected.csv")
