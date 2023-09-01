library(depthtestr)
library(tidyverse)
library(patchwork)
library(colorBlindness)

# Overview: this script runs the 2-d uniform ball example: Figure 2 
# in section 3.4 of the manuscript (Evans and Berenhaut)

# Helper functions used in the depthtestr package to calculate 
# local community depth (LCD) and pairwise depth (PD)
# Also need to source the C++ file
Rcpp::sourceCpp("depth_helpers.cpp")

lcdDepthObs <- function(sampB_obs, sampA, dists){
  N1 = nrow(sampA)
  
  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))
  
  return(1/N1 * lcdDepthDists(new_dists, dists, N1))
}

## compute the local community depth of points in sampB,
## relative to sampA (the reference sample)
## Returns a vector of the same length as sampB, with a 
## depth for each point in sampB
depthFunLCD <- function(sampB, sampA){
  dists <- as.matrix(stats::dist(sampA))
  apply(sampB, 1, function(x) lcdDepthObs(x, sampA, dists))
}


pdDepthObs <- function(sampB_obs, sampA, dists){
  N1 = nrow(sampA)
  
  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))
  
  return(1/(N1 * (N1 - 1)) * pdDepthDists(new_dists, dists, N1))
}

## compute the pairwise depth of points in sampB,
## relative to sampA (the reference sample)
## Returns a vector of the same length as sampB, with a 
## depth for each point in sampB
depthFunPD <- function(sampB, sampA){
  dists <- as.matrix(stats::dist(sampA))
  apply(sampB, 1, function(x) pdDepthObs(x, sampA, dists))
}


# Helper function for generating data from a uniform ball
# n = sample size
# d = dimension
# c = center of ball
# r = radius of ball
runif_ball <- function(n, d, c, r){
  radii <- runif(n)^(1/d)
  z = matrix(rnorm(n*d),nrow=n,ncol=d)
  x = r*z/rep(sqrt(rowSums(z^2))/radii,1)
  sweep(x, 2, c, "+")
}



# Part 1: Plot local community depth and pairwise depth 
# for a mixture of three uniform balls, for different amounts of
# separation between the balls


set.seed(712)

N1 = 1000 # size of reference sample used to estimate depth
d <- 2 # dimension of data
alpha = c(1/3, 1/3, 1/3) # mixing proportions for the three balls


## begin with distance 1 between any two centers
sep = 1

## reference sample for estimating depth
y1 <- sample(c(0, 1, 2), N1, replace = T, prob = alpha)
samp1 <- runif_ball(N1, d = d, c = c(0, 0), r = 0.5)*(y1 == 0) +
  runif_ball(N1, d = d, c = c(0.5*sep, 0.5*sep*sqrt(3)), r = 0.5)*(y1 == 1) +
  runif_ball(N1, d = d, c = c(sep, 0), r = 0.5)*(y1 == 2)

## grid of points at which to evaluate depth
samp2 = expand.grid(seq(-2, 3, 0.05),
                    seq(-2, 3, 0.05))

## evaluate LCD at each point on the grid
lcd_depths <- depthFunLCD(samp2, samp1)

p1 <- data.frame(x1 = samp2[,1], x2 = samp2[,2], 
           depth = lcd_depths) |> 
  ggplot(aes(x = x1, y = x2, fill = depth)) +
  geom_raster() +
  scale_fill_gradientn(colors=rev(ModifiedSpectralScheme11Steps),
                       limits=c(0, 0.81)) +
  theme_classic() +
  labs(x = "", y = "", fill = "Depth", 
       title = "Local community depth, separation = 1")

## evaluate PD at each point on the grid
pd_depths <- depthFunPD(samp2, samp1)

p2 <- data.frame(x1 = samp2[,1], x2 = samp2[,2], 
           depth = pd_depths) |> 
  ggplot(aes(x = x1, y = x2, fill = depth)) +
  geom_raster() +
  scale_fill_gradientn(colors=rev(ModifiedSpectralScheme11Steps),
                       limits=c(0, 0.81)) +
  theme_classic() +
  labs(x = "", y = "", fill = "Depth", 
       title = "Pairwise depth, separation = 1")




## Repeat again for distance 2 between centers
sep = 2

## reference sample for estimating depth
y1 <- sample(c(0, 1, 2), N1, replace = T, prob = alpha)
samp1 <- runif_ball(N1, d = d, c = c(0, 0), r = 0.5)*(y1 == 0) +
  runif_ball(N1, d = d, c = c(0.5*sep, 0.5*sep*sqrt(3)), r = 0.5)*(y1 == 1) +
  runif_ball(N1, d = d, c = c(sep, 0), r = 0.5)*(y1 == 2)

## grid of points at which to evaluate depth
samp2 = expand.grid(seq(-3, 5, 0.1),
                    seq(-3, 5, 0.1))

## evaluate LCD at each point on the grid
lcd_depths <- depthFunLCD(samp2, samp1)

p3 <- data.frame(x1 = samp2[,1], x2 = samp2[,2], 
                 depth = lcd_depths) |> 
  ggplot(aes(x = x1, y = x2, fill = depth)) +
  geom_raster() +
  scale_fill_gradientn(colors=rev(ModifiedSpectralScheme11Steps),
                       limits=c(0, 0.81)) +
  theme_classic() +
  labs(x = "", y = "", fill = "Depth", 
       title = "Local community depth, separation = 2")

## evaluate PD at each point on the grid
pd_depths <- depthFunPD(samp2, samp1)

p4 <- data.frame(x1 = samp2[,1], x2 = samp2[,2], 
                 depth = pd_depths) |> 
  ggplot(aes(x = x1, y = x2, fill = depth)) +
  geom_raster() +
  scale_fill_gradientn(colors=rev(ModifiedSpectralScheme11Steps),
                       limits = c(0, 0.81)) +
  theme_classic() +
  labs(x = "", y = "", fill = "Depth", 
       title = "Pairwise depth, separation = 2")

## repeat again for distance 3 between centers
sep = 3

y1 <- sample(c(0, 1, 2), N1, replace = T, prob = alpha)
samp1 <- runif_ball(N1, d = d, c = c(0, 0), r = 0.5)*(y1 == 0) +
  runif_ball(N1, d = d, c = c(0.5*sep, 0.5*sep*sqrt(3)), r = 0.5)*(y1 == 1) +
  runif_ball(N1, d = d, c = c(sep, 0), r = 0.5)*(y1 == 2)

samp2 = expand.grid(seq(-4, 7, 0.1),
                    seq(-4, 7, 0.1))

lcd_depths <- depthFunLCD(samp2, samp1)

p5 <- data.frame(x1 = samp2[,1], x2 = samp2[,2], 
                 depth = lcd_depths) |> 
  ggplot(aes(x = x1, y = x2, fill = depth)) +
  geom_raster() +
  scale_fill_gradientn(colors=rev(ModifiedSpectralScheme11Steps),
                       limits = c(0, 0.81)) +
  theme_classic() +
  labs(x = "", y = "", fill = "Depth", 
       title = "Local community depth, separation = 3")


pd_depths <- depthFunPD(samp2, samp1)

p6 <- data.frame(x1 = samp2[,1], x2 = samp2[,2], 
                 depth = pd_depths) |> 
  ggplot(aes(x = x1, y = x2, fill = depth)) +
  geom_raster() +
  scale_fill_gradientn(colors=rev(ModifiedSpectralScheme11Steps),
                       limits = c(0, 0.81)) +
  theme_classic() +
  labs(x = "", y = "", fill = "Depth", 
       title = "Pairwise depth, separation = 3")




# Part 2: Calculating power for different distances between the uniform balls
# Data X_1,...,X_N1 and Y1,...,Y_N2 are sampled from distributions F and G,
# respectively. Here F and G are both mixtures of three uniform balls
# The radii in F are all 0.5
# The radii in G are 0.5 + delta/sqrt(N1 + N2)
# We calculate power for testing F = G

set.seed(279)

N1 = 100 # size of first sample
N2 = 100 # size of second sample
nsamp <- 500 # number of repetitions to estimate power

delta <- 2 
separations <- c(1, 1.5, 2, 2.5, 3) # distances between centers

# vectors to store the results
power_lcd <- c()
power_pd <- c()

for(j in 1:length(separations)){
  
  stats_lcd <- c()
  stats_null_lcd <- c()
  stats_pd <- c()
  stats_null_pd <- c()
  
  pvals_lcd <- c()
  pvals_pd <- c()
  
  sep <- separations[j]
  
  for(i in 1:nsamp){
    
    ## Sample from F
    y1 <- sample(c(0, 1, 2), N1, replace = T, prob = alpha)
    samp1 <- runif_ball(N1, d = d, c = c(0, 0), r = 0.5)*(y1 == 0) +
      runif_ball(N1, d = d, c = c(0.5*sep, 0.5*sep*sqrt(3)), r = 0.5)*(y1 == 1) +
      runif_ball(N1, d = d, c = c(sep, 0), r = 0.5)*(y1 == 2)
    
    
    ## Sample from G
    y2 <- sample(c(0, 1, 2), N2, replace = T, prob = alpha)
    samp2 <- runif_ball(N2, d = d, c = c(0, 0), 
                        r = 0.5 + delta/sqrt(N1 + N2))*(y2 == 0) +
      runif_ball(N2, d = d, c = c(0.5*sep, 0.5*sep*sqrt(3)), 
                 r = 0.5 + delta/sqrt(N1 + N2))*(y2 == 1) +
      runif_ball(N2, d = d, c = c(sep, 0), 
                 r = 0.5 + delta/sqrt(N1 + N2))*(y2 == 2)
    
    ## A second sample from F, used to calculate power while 
    ## making type I error comparable (approximates empirical distribution
    ## of the p-value)
    y3 <- sample(c(0, 1, 2), N1, replace = T, prob = alpha)
    samp3 <- runif_ball(N1, d = d, c = c(0, 0), r = 0.5)*(y3 == 0) +
      runif_ball(N1, d = d, c = c(0.5*sep, 0.5*sep*sqrt(3)), r = 0.5)*(y3 == 1) +
      runif_ball(N1, d = d, c = c(sep, 0), r = 0.5)*(y3 == 2)
    
    test_lcd_12 <- suppressMessages(depthTest(samp1, samp2, "lcd"))
    test_pd_12 <- suppressMessages(depthTest(samp1, samp2, "pd"))
    
    test_lcd_13 <- suppressMessages(depthTest(samp1, samp3, "lcd"))
    test_pd_13 <- suppressMessages(depthTest(samp1, samp3, "pd"))
    
    stats_lcd[i] <- test_lcd_12$teststat
    stats_null_lcd[i] <- test_lcd_13$teststat
    stats_pd[i] <- test_pd_12$teststat
    stats_null_pd[i] <- test_pd_13$teststat
    
    pvals_lcd[i] <- test_lcd_12$pval
    pvals_pd[i] <- test_pd_12$pval
    
    print(i)
  }
  
  # power of LCD
  power_lcd[j] <- mean(stats_lcd > quantile(stats_null_lcd, 0.975)) +
    mean(stats_lcd < quantile(stats_null_lcd, 0.025))
  
  # power of PD
  power_pd[j] <- mean(stats_pd > quantile(stats_null_pd, 0.975)) +
    mean(stats_pd < quantile(stats_null_pd, 0.025))
}

## plot power against separation
p7 <- data.frame(sep = rep(separations, 2),
           method = rep(c("LCD", "PD"), each=length(separations)),
           power = c(power_lcd, power_pd)) |>
  ggplot(aes(x = sep, y = power, color = method, shape = method)) +
  geom_point(size=2) +
  geom_line(lwd=0.8) +
  theme_bw() +
  labs(x = "Separation", y = "Power",
       color = "Depth", shape = "Depth")
  

# Finally, combine all the plots into a single figure

layout <- "
AABBCC##
AABBCCDD
EEFFGGDD
EEFFGG##
"

pdf(file = "example_2_plot.pdf",
    width = 16, height = 8)

p1 + p3 + p5 + p7 + p2 + p4 + p6 + 
  plot_layout(design = layout,
              guides = "collect") & theme(legend.position = 'bottom')

dev.off()

