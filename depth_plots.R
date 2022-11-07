library(fda.usc)
library(tidyverse)
library(gridExtra)

weight_fun <- function(x1, x2, cdf){
  upper <- max(2*x1 - x2, 2*x2 - x1)
  lower <- min(2*x1 - x2, 2*x2 - x1)
  return(1 - cdf(upper) + cdf(lower))
}

cdf_fun <- function(x){
  pnorm(x)*alpha + pnorm(x, mean=mu2, sd=sd2)*(1 - alpha)
}

alpha = 0.5
nsamp = 20000
mu2 = 5
sd2 = 0.25
labels1 <- rbinom(nsamp, 1, alpha)
labels2 <- rbinom(nsamp, 1, alpha)
samp1 <- rnorm(nsamp)*labels1 + rnorm(nsamp, mean=mu2, sd=sd2)*(1 - labels1)
samp2 <- rnorm(nsamp)*labels2 + rnorm(nsamp, mean=mu2, sd=sd2)*(1 - labels2)

weights <- c()
for(i in 1:length(samp1)){
  weights[i] <- weight_fun(samp1[i], samp2[i], cdf_fun)
}

samp_dists <- abs(samp1 - samp2)

depth_fun <- function(y){
  y_dists_1 <- abs(y - samp1)
  y_dists_2 <- abs(y - samp2)
  mean(weights * (y_dists_1 < samp_dists | y_dists_2 < samp_dists))
}

points <- seq(-5, 9, 0.02)

labels <- rbinom(nsamp, 1, alpha)
samp_points <- rnorm(nsamp)*labels + rnorm(nsamp, mean=mu2, sd=sd2)*(1 - labels)


depths <- sapply(points, depth_fun)
depths_samp <- sapply(samp_points, depth_fun)
r_y <- sapply(depths, function(y){mean(depths_samp <= y)})


depths_dens <- dnorm(points)*alpha + dnorm(points, mu2, sd2)*(1-alpha)
depths_samp_dens <- dnorm(samp_points)*alpha + dnorm(samp_points, mu2, sd2)*(1-alpha)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})

deriv <- -1*(alpha * dnorm(points, 0, 1) + 
  (1 - alpha) * dnorm(points, mu2, sd2) - 
  alpha*(points - 0)^2 * dnorm(points, 0, 1) - 
  (1-alpha)*(points - mu2)^2/(sd2^2) * dnorm(points, mu2, sd2))

p1 <- data.frame(x = points, y = depths_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Density", title = "Gaussian mixture") +
  theme(text = element_text(size = 15))

p2 <- data.frame(x = points, y = depths) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p3 <- data.frame(x = points, y = md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p4 <- data.frame(x = points, y = deriv) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Derivative", title = "Derivative wrt scale parameter") +
  theme(text = element_text(size = 15))

p5 <- data.frame(x = points, y = r_y) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p6 <- data.frame(x = points, y = r_y_md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

pdf("mixture_density_example.pdf", width=15, height=6)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)
dev.off()

mean(-1*deriv * r_y_md)
mean(-1*deriv * r_y)
mean(-1*deriv * r_y_dens)

