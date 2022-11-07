library(tidyverse)
library(gridExtra)

power_vs_distance <- read.csv("mixture_results_power_vs_distance.csv")

p1 <- data.frame(power = c(power_vs_distance$power_adjusted_mah,
                           power_vs_distance$power_adjusted_hd,
                           power_vs_distance$power_adjusted_pd,
                           power_vs_distance$power_adjusted_lcd,
                           power_vs_distance$power_adjusted_lp,
                           power_vs_distance$power_adjusted_pvb,
                           power_vs_distance$power_cf),
                 method = rep(c("MD", "HD", "PD", "LCD", "Lp", "PVB", "CF"),
                              each = 8),
                 distance = rep(c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), 7)) %>%
  ggplot(aes(x = distance, y = power, color = method, shape=method)) +
  geom_point(size=2) +
  geom_line(lwd=0.75) +
  theme_bw() +
  labs(x = "Distance between means",
       y = "Power",
       color = "Method", shape = "Method") +
  theme(text = element_text(size = 15))


power_vs_sd <- read.csv("mixture_results_power_vs_sd.csv")

p2 <- data.frame(power = c(power_vs_sd$power_adjusted_mah,
                           power_vs_sd$power_adjusted_hd,
                           power_vs_sd$power_adjusted_pd,
                           power_vs_sd$power_adjusted_lcd,
                           power_vs_sd$power_adjusted_lp,
                           power_vs_sd$power_adjusted_pvb,
                           power_vs_sd$power_cf),
                 method = rep(c("MD", "HD", "PD", "LCD", "Lp", "PVB", "CF"),
                              each = 7),
                 sd = seq(0.25, 1, 0.125)) %>%
  ggplot(aes(x = sd, y = power, color = method, shape=method)) +
  geom_point(size=2) +
  geom_line(lwd=0.75) +
  theme_bw() +
  labs(x = "Ratio of standard deviations",
       y = "Power",
       color = "Method", shape = "Method") +
  theme(text = element_text(size = 15))

pdf("mixture_performance_power.pdf", width = 10, height=3)
grid.arrange(p1, p2, ncol=2)
dev.off()
