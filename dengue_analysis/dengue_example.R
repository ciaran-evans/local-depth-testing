library(tidyverse)
library(gTests)
library(depthtestr)
library(vegan)
library(mvtnorm)


dengue_original <- read.csv("dengue_original.csv")


dengue <- dengue_original %>%
  select(SiteNo, DayDisease, Temp, HCT, PLT, 
         NEUCount, LYMCount, MONOCount,
         ALB, AST, ALT, CK,
         Lab_Confirmed_Dengue) %>%
  rename(SiteNumber = SiteNo, 
         DiseaseDay = DayDisease,
         Temperature = Temp, 
         Dengue = Lab_Confirmed_Dengue) %>%
  mutate(Dengue = 2 - Dengue,
         SiteNumber = ifelse(SiteNumber == 15, 5, 
                             ifelse(SiteNumber == 16, 6, SiteNumber))) %>% 
  drop_na()





nreps <- 500
n_dengue <- seq(0, 50, 5)

power_hd <- c()
power_mah <- c()
power_lp <- c()
power_pd <- c()
power_lcd <- c()
power_cf <- c()
power_pvb <- c()

power_adjusted_hd <- c()
power_adjusted_mah <- c()
power_adjusted_pvb <- c()
power_adjusted_lp <- c()
power_adjusted_lcd <- c()
power_adjusted_pd <- c()


dengue_0 <- dengue %>%
  mutate(Temperature = scale(Temperature),
         PLT = scale(PLT),
         HCT = scale(HCT),
         ALB = scale(ALB),
         AST = scale(AST),
         CK = scale(CK),
         NEUCount = scale(NEUCount),
         LYMCount = scale(LYMCount),
         MONOCount = scale(MONOCount)) %>%
  filter(Dengue == 0) %>%
  select(Temperature, PLT, HCT, 
         ALB, AST, CK,
         NEUCount, LYMCount, MONOCount) 

dengue_1 <- dengue %>%
  mutate(Temperature = scale(Temperature),
         PLT = scale(PLT),
         HCT = scale(HCT),
         ALB = scale(ALB),
         AST = scale(AST),
         CK = scale(CK),
         NEUCount = scale(NEUCount),
         LYMCount = scale(LYMCount),
         MONOCount = scale(MONOCount)) %>%
  filter(Dengue == 1) %>%
  select(Temperature, PLT, HCT, 
         ALB, AST, CK,
         NEUCount, LYMCount, MONOCount) 

for(j in 1:length(n_dengue)){
  pvals_hd <- c()
  pvals_mah <- c()
  pvals_lp <- c()
  pvals_pd <- c()
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
  stats_lcd <- c()
  stats_null_lcd <- c()
  stats_pd <- c()
  stats_null_pd <- c()
  
  num_positive <- n_dengue[j]
  
  set.seed(j)
  
  for(i in 1:nreps){
    
    samp1 <- dengue_0 %>%
      slice_sample(n = 50)
    
    samp2 <- dengue_1 %>%
      slice_sample(n = num_positive) %>%
      rbind(dengue_0 %>%
              slice_sample(n = (50 - num_positive)))
    
    samp3 <- dengue_0 %>%
      slice_sample(n = 50)
    
    # test_hd_12 <- depthTest(samp1, samp2, "halfspace")
    # test_mah_12 <- depthTest(samp1, samp2, "mahalanobis")
    # test_pvb_12 <- depthTest(samp1, samp2, "pvb")
    # test_lp_12 <- depthTest(samp1, samp2, "lp")
    # test_lcd_12 <- depthTest(samp1, samp2, "lcd")
    # test_pd_12 <- depthTest(samp1, samp2, "pd")
    # 
    # test_hd_13 <- depthTest(samp1, samp3, "halfspace")
    # test_mah_13 <- depthTest(samp1, samp3, "mahalanobis")
    # test_pvb_13 <- depthTest(samp1, samp3, "pvb")
    # test_lp_13 <- depthTest(samp1, samp3, "lp")
    # test_lcd_13 <- depthTest(samp1, samp3, "lcd")
    # test_pd_13 <- depthTest(samp1, samp3, "pd")
    
    test_hd_12 <- suppressMessages(depthTest(samp2, samp1, "halfspace"))
    test_mah_12 <- suppressMessages(depthTest(samp2, samp1, "mahalanobis"))
    test_pvb_12 <- suppressMessages(depthTest(samp2, samp1, "pvb"))
    test_lp_12 <- suppressMessages(depthTest(samp2, samp1, "lp"))
    test_lcd_12 <- suppressMessages(depthTest(samp2, samp1, "lcd"))
    test_pd_12 <- suppressMessages(depthTest(samp2, samp1, "pd"))
    
    test_hd_13 <- suppressMessages(depthTest(samp3, samp1, "halfspace"))
    test_mah_13 <- suppressMessages(depthTest(samp3, samp1, "mahalanobis"))
    test_pvb_13 <- suppressMessages(depthTest(samp3, samp1, "pvb"))
    test_lp_13 <- suppressMessages(depthTest(samp3, samp1, "lp"))
    test_lcd_13 <- suppressMessages(depthTest(samp3, samp1, "lcd"))
    test_pd_13 <- suppressMessages(depthTest(samp3, samp1, "pd"))
    
    stats_hd[i] <- test_hd_12$teststat
    stats_null_hd[i] <- test_hd_13$teststat
    stats_mah[i] <- test_mah_12$teststat
    stats_null_mah[i] <- test_mah_13$teststat
    stats_pvb[i] <- test_pvb_12$teststat
    stats_null_pvb[i] <- test_pvb_13$teststat
    stats_lp[i] <- test_lp_12$teststat
    stats_null_lp[i] <- test_lp_13$teststat
    stats_lcd[i] <- test_lcd_12$teststat
    stats_null_lcd[i] <- test_lcd_13$teststat
    stats_pd[i] <- test_pd_12$teststat
    stats_null_pd[i] <- test_pd_13$teststat
    
    pvals_hd[i] <- test_hd_12$pval
    pvals_mah[i] <- test_mah_12$pval
    pvals_lp[i] <- test_lp_12$pval
    pvals_pd[i] <- test_pd_12$pval
    pvals_lcd[i] <- test_lcd_12$pval
    pvals_pvb[i] <- test_pvb_12$pval
    
    N1 <- 50
    N2 <- 50
    
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
  power_adjusted_lcd[j] <- mean(stats_lcd > quantile(stats_null_lcd, 0.975)) + 
    mean(stats_lcd < quantile(stats_null_lcd, 0.025))
  
  print(num_positive)
}

output <- data.frame(num_positive = n_dengue,
                     power_hd = power_hd,
                     power_mah = power_mah,
                     power_lp = power_lp,
                     power_pd = power_pd,
                     power_lcd = power_lcd,
                     power_pvb = power_pvb,
                     power_cf = power_cf,
                     power_adjusted_hd = power_adjusted_hd,
                     power_adjusted_mah = power_adjusted_mah,
                     power_adjusted_lp = power_adjusted_lp,
                     power_adjusted_pd = power_adjusted_pd,
                     power_adjusted_lcd = power_adjusted_lcd,
                     power_adjusted_pvb = power_adjusted_pvb)

write_csv(output, file = "dengue_power_vs_positives.csv")



# using a subset of features with the biggest scale change

nreps <- 500
n_dengue <- seq(0, 50, 5)

power_hd <- c()
power_mah <- c()
power_lp <- c()
power_pd <- c()
power_lcd <- c()
power_cf <- c()
power_pvb <- c()

power_adjusted_hd <- c()
power_adjusted_mah <- c()
power_adjusted_pvb <- c()
power_adjusted_lp <- c()
power_adjusted_lcd <- c()
power_adjusted_pd <- c()

dengue_0 <- dengue %>%
  mutate(Temperature = scale(Temperature),
         PLT = scale(PLT),
         HCT = scale(HCT),
         ALB = scale(ALB),
         AST = scale(AST),
         CK = scale(CK),
         NEUCount = scale(NEUCount),
         LYMCount = scale(LYMCount),
         MONOCount = scale(MONOCount)) %>%
  filter(Dengue == 0) %>%
  select(NEUCount, LYMCount, MONOCount) 

dengue_1 <- dengue %>%
  mutate(Temperature = scale(Temperature),
         PLT = scale(PLT),
         HCT = scale(HCT),
         ALB = scale(ALB),
         AST = scale(AST),
         CK = scale(CK),
         NEUCount = scale(NEUCount),
         LYMCount = scale(LYMCount),
         MONOCount = scale(MONOCount)) %>%
  filter(Dengue == 1) %>%
  select(NEUCount, LYMCount, MONOCount) 

for(j in 1:length(n_dengue)){
  pvals_hd <- c()
  pvals_mah <- c()
  pvals_lp <- c()
  pvals_pd <- c()
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
  stats_lcd <- c()
  stats_null_lcd <- c()
  stats_pd <- c()
  stats_null_pd <- c()
  
  num_positive <- n_dengue[j]
  
  set.seed(j)
  
  for(i in 1:nreps){
    
    samp1 <- dengue_0 %>%
      slice_sample(n = 50)
    
    samp2 <- dengue_1 %>%
      slice_sample(n = num_positive) %>%
      rbind(dengue_0 %>%
              slice_sample(n = (50 - num_positive)))
    
    samp3 <- dengue_0 %>%
      slice_sample(n = 50)
    
    # test_hd_12 <- depthTest(samp1, samp2, "halfspace")
    # test_mah_12 <- depthTest(samp1, samp2, "mahalanobis")
    # test_pvb_12 <- depthTest(samp1, samp2, "pvb")
    # test_lp_12 <- depthTest(samp1, samp2, "lp")
    # test_lcd_12 <- depthTest(samp1, samp2, "lcd")
    # test_pd_12 <- depthTest(samp1, samp2, "pd")
    # 
    # test_hd_13 <- depthTest(samp1, samp3, "halfspace")
    # test_mah_13 <- depthTest(samp1, samp3, "mahalanobis")
    # test_pvb_13 <- depthTest(samp1, samp3, "pvb")
    # test_lp_13 <- depthTest(samp1, samp3, "lp")
    # test_lcd_13 <- depthTest(samp1, samp3, "lcd")
    # test_pd_13 <- depthTest(samp1, samp3, "pd")
    
    test_hd_12 <- suppressMessages(depthTest(samp2, samp1, "halfspace"))
    test_mah_12 <- suppressMessages(depthTest(samp2, samp1, "mahalanobis"))
    test_pvb_12 <- suppressMessages(depthTest(samp2, samp1, "pvb"))
    test_lp_12 <- suppressMessages(depthTest(samp2, samp1, "lp"))
    test_lcd_12 <- suppressMessages(depthTest(samp2, samp1, "lcd"))
    test_pd_12 <- suppressMessages(depthTest(samp2, samp1, "pd"))
    
    test_hd_13 <- suppressMessages(depthTest(samp3, samp1, "halfspace"))
    test_mah_13 <- suppressMessages(depthTest(samp3, samp1, "mahalanobis"))
    test_pvb_13 <- suppressMessages(depthTest(samp3, samp1, "pvb"))
    test_lp_13 <- suppressMessages(depthTest(samp3, samp1, "lp"))
    test_lcd_13 <- suppressMessages(depthTest(samp3, samp1, "lcd"))
    test_pd_13 <- suppressMessages(depthTest(samp3, samp1, "pd"))
    
    stats_hd[i] <- test_hd_12$teststat
    stats_null_hd[i] <- test_hd_13$teststat
    stats_mah[i] <- test_mah_12$teststat
    stats_null_mah[i] <- test_mah_13$teststat
    stats_pvb[i] <- test_pvb_12$teststat
    stats_null_pvb[i] <- test_pvb_13$teststat
    stats_lp[i] <- test_lp_12$teststat
    stats_null_lp[i] <- test_lp_13$teststat
    stats_lcd[i] <- test_lcd_12$teststat
    stats_null_lcd[i] <- test_lcd_13$teststat
    stats_pd[i] <- test_pd_12$teststat
    stats_null_pd[i] <- test_pd_13$teststat
    
    pvals_hd[i] <- test_hd_12$pval
    pvals_mah[i] <- test_mah_12$pval
    pvals_lp[i] <- test_lp_12$pval
    pvals_pd[i] <- test_pd_12$pval
    pvals_lcd[i] <- test_lcd_12$pval
    pvals_pvb[i] <- test_pvb_12$pval
    
    N1 <- 50
    N2 <- 50
    
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
  power_adjusted_lcd[j] <- mean(stats_lcd > quantile(stats_null_lcd, 0.975)) + 
    mean(stats_lcd < quantile(stats_null_lcd, 0.025))
  
  print(num_positive)
}

output <- data.frame(num_positive = n_dengue,
                     power_hd = power_hd,
                     power_mah = power_mah,
                     power_lp = power_lp,
                     power_pd = power_pd,
                     power_lcd = power_lcd,
                     power_pvb = power_pvb,
                     power_cf = power_cf,
                     power_adjusted_hd = power_adjusted_hd,
                     power_adjusted_mah = power_adjusted_mah,
                     power_adjusted_lp = power_adjusted_lp,
                     power_adjusted_pd = power_adjusted_pd,
                     power_adjusted_lcd = power_adjusted_lcd,
                     power_adjusted_pvb = power_adjusted_pvb)

write_csv(output, file = "dengue_power_vs_positives_subset.csv")



# plot results

dengue_power_vs_positives <- read_csv("dengue_power_vs_positives.csv")
dengue_power_vs_positives_subset <- read_csv("dengue_power_vs_positives_subset.csv")

p1 <- dengue %>%
  mutate(Dengue = ifelse(Dengue == 0, "negative", "positive")) %>%
  ggplot(aes(x = LYMCount, y = NEUCount, color = Dengue)) +
  geom_density_2d(contour_var = "ndensity", lwd=0.7) +
  theme_bw() +
  labs(x = "Lymphocyte count", 
       y = "Neutrophil count",
       color = "Dengue status") +
  theme(text = element_text(size = 15))

p2 <- dengue %>%
  mutate(Dengue = ifelse(Dengue == 0, "negative", "positive")) %>%
  ggplot(aes(x = Temperature, y = HCT, color = Dengue)) +
  geom_density_2d(contour_var = "ndensity", lwd=0.7) +
  theme_bw() +
  labs(x = "Temperature", 
       y = "Hematocrit",
       color = "Dengue status") +
  theme(text = element_text(size = 15))


p3 <- data.frame(n = rep(dengue_power_vs_positives$num_positive, 7),
           power = c(dengue_power_vs_positives$power_adjusted_hd,
                     dengue_power_vs_positives$power_adjusted_mah,
                     dengue_power_vs_positives$power_adjusted_lp,
                     dengue_power_vs_positives$power_adjusted_pd,
                     dengue_power_vs_positives$power_adjusted_lcd,
                     dengue_power_vs_positives$power_adjusted_pvb,
                     dengue_power_vs_positives$power_cf),
           type = rep(c("HD", "MD", "Lp", "PD", "LCD", "PVB", "CF"), each=11)) %>%
  ggplot(aes(x = n, y = power, color = type, shape=type)) +
  geom_point(size=2) +
  geom_line(lwd=0.75) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  labs(x = "Number of positive cases",
       y = "Power",
       title = "Test using all quantitative features",
       color = "Method", shape = "Method") +
  scale_y_continuous(limits=c(0, 1))

p4 <- data.frame(n = rep(dengue_power_vs_positives_subset$num_positive, 7),
                 power = c(dengue_power_vs_positives_subset$power_adjusted_hd,
                           dengue_power_vs_positives_subset$power_adjusted_mah,
                           dengue_power_vs_positives_subset$power_adjusted_lp,
                           dengue_power_vs_positives_subset$power_adjusted_pd,
                           dengue_power_vs_positives_subset$power_adjusted_lcd,
                           dengue_power_vs_positives_subset$power_adjusted_pvb,
                           dengue_power_vs_positives_subset$power_cf),
                 type = rep(c("HD", "MD", "Lp", "PD", "LCD", "PVB", "CF"), each=11)) %>%
  ggplot(aes(x = n, y = power, color = type, shape=type)) +
  geom_point(size=2) +
  geom_line(lwd=0.75) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  labs(x = "Number of positive cases",
       y = "Power",
       title = "Test using a subset of features",
       color = "Method", shape = "Method") +
  scale_y_continuous(limits=c(0,1))

pdf("dengue_features.pdf", width = 10, height = 3)
grid.arrange(p1, p2, ncol=2)
dev.off()

pdf("dengue_results.pdf", width=10, height=3)
grid.arrange(p3, p4, ncol=2)
dev.off()


