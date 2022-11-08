library(tidyverse)
library(gridExtra)

results_files <- c("normal_location_shift_results.csv",
                   "normal_scale_shift_results.csv",
                   "lognormal_shift_results.csv",
                   "mixture_results_2_components.csv",
                   "mixture_results_3_components.csv")

lcd_error <- c()
lcd_power <- c()
pald_error <- c()
pald_power <- c()
scenario <- c(rep("Normal location shift", 5),
              rep("Normal scale shift", 5),
              rep("Lognormal shift", 5),
              rep("Mixture, 2 components", 5),
              rep("Mixture, 3 components", 5))

for(i in 1:5){
  res <- read.csv(results_files[i])
  lcd_error <- c(lcd_error, res$type_1_error_lcd)
  lcd_power <- c(lcd_power, res$power_adjusted_lcd)
  pald_error <- c(pald_error, res$type_1_error_pald)
  pald_power <- c(pald_power, res$power_adjusted_pald)
}

p1 <- data.frame(lcd_error, lcd_power, pald_error, pald_power, scenario) %>%
  ggplot(aes(x = lcd_error, y = pald_error, color = scenario,
             shape = scenario)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(x = "LCD type I error", y = "PaLD type I error",
       color = "Scenario", shape = "Scenario")

p2 <- data.frame(lcd_error, lcd_power, pald_error, pald_power, scenario) %>%
  ggplot(aes(x = lcd_power, y = pald_power, color = scenario,
             shape = scenario)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(x = "LCD power", y = "PaLD power",
       color = "Scenario", shape = "Scenario")

pdf("lcd_vs_pald_results.pdf", 10, 4)
grid.arrange(p1, p2, ncol=2)
dev.off()