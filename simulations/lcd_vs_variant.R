library(tidyverse)
library(patchwork)

results_files <- c("normal_scale_shift_results.csv",
                   #"normal_location_shift_results.csv",
                   "lognormal_shift_results.csv",
                   "mixture_results_2_components.csv",
                   "mixture_results_3_components.csv",
                   "uniform_ball_scale_shift_results.csv",
                   "uniform_ball_scale_shift_3_components_results.csv")

lcd_error <- c()
lcd_power <- c()
lcd_variant_error <- c()
lcd_variant_power <- c()
scenario <- c(rep("Normal scale shift", 5),
              #rep("Normal location shift", 5),
              rep("Lognormal shift", 5),
              rep("Gaussian mixture, 2 components", 5),
              rep("Gaussian mixture, 3 components", 5),
              rep("Uniform ball", 5),
              rep("Uniform ball mixture, 3 components", 5))

for(i in 1:length(results_files)){
  res <- read.csv(results_files[i])
  lcd_error <- c(lcd_error, res$type_1_error_lcd)
  lcd_power <- c(lcd_power, res$power_adjusted_lcd)
  lcd_variant_error <- c(lcd_variant_error, res$type_1_error_lcd_variant)
  lcd_variant_power <- c(lcd_variant_power, res$power_adjusted_lcd_variant)
}

p1 <- data.frame(lcd_error, lcd_power, lcd_variant_error, 
                 lcd_variant_power, scenario) %>%
  ggplot(aes(x = lcd_error, y = lcd_variant_error, color = scenario,
             shape = scenario)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(x = "LCD type I error", y = "LCD variant type I error",
       color = "Scenario", shape = "Scenario")

p2 <- data.frame(lcd_error, lcd_power, lcd_variant_error, 
                 lcd_variant_power, scenario) %>%
  ggplot(aes(x = lcd_power, y = lcd_variant_power, color = scenario,
             shape = scenario)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(x = "LCD power", y = "LCD variant power",
       color = "Scenario", shape = "Scenario")

pdf("lcd_vs_variant_results.pdf", 10, 4)
combined <- p1 + p2 & theme(legend.position = "right")
combined + plot_layout(guides = "collect")
dev.off()
