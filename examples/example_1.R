
cdf_fun <- function(x){
  if(length(x) > 1){
    rowSums(sapply(1:length(alpha), function(i) alpha[i]*pnorm(x, mu[i], sd[i])))
  } else {
    sum(sapply(1:length(alpha), function(i) alpha[i]*pnorm(x, mu[i], sd[i])))
  }
}

density_fun <- function(x){
  rowSums(sapply(1:length(alpha), function(i) alpha[i]*dnorm(x, mu[i], sd[i])))
}

deriv_fun <- function(x){
  -1*(rowSums(sapply(1:length(alpha), function(i) alpha[i]*dnorm(x, mu[i], sd[i]))) -
    rowSums(sapply(1:length(alpha), 
                   function(i) alpha[i]*(x - mu[i])^2/(sd[i]^2)*dnorm(x, mu[i], sd[i]))))
}

deriv_fun_2 <- function(x){
  dnorm(x, mu[1], sd[1]) - dnorm(x, mu[2], sd[2])
}

deriv_fun_3 <- function(x){
  rowSums(sapply(1:length(alpha), 
                 function(i) alpha[i] * dnorm(x, mu[i], sd[i]) * (x - mu[i])/(sd[i]^2)))
}

data_fun <- function(nsamp){
  labels <- sample(1:length(alpha), nsamp, replace=T, prob = alpha)
  rnorm(nsamp, mu[labels], sd[labels])
}

weight_fun <- function(y, x, cdf){
  if(y > x){
    (cdf(2*y - x) - cdf((x+y)/2))/(cdf(2*y - x) - cdf(2*x - y))
  } else {
    (cdf((x+y)/2) - cdf(2*y - x))/(cdf(2*x - y) - cdf(2*y - x))
  }
}

local_depth_fun <- function(y){
  weights <- sapply(1:length(samp1[1:1000]), 
                    function(i) weight_fun(y, samp1[1:1000][i], cdf_fun))
  mean(weights)
}



set.seed(3)

mu <- c(0, 0)
alpha <- c(0.5, 0.5)
sd <- c(1, 1)


points <- seq(-5, 5, 0.02)
derivs <- deriv_fun(points)
derivs_3 <- deriv_fun_3(points)


nsamp = 20000
samp1 <- data_fun(nsamp)
samp2 <- data_fun(nsamp)
samp_points <- data_fun(nsamp)

depths_dens <- density_fun(points)
depths_samp_dens <- density_fun(samp_points)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})


local_depths <- rep(NA, length(points))
for(i in 1:length(points)){
  local_depths[i] <- local_depth_fun(points[i])
  print(i)
}
lcd_smooth <- loess(local_depths ~ points, span = 0.05)
local_depths_samp <- predict(lcd_smooth, samp_points)
r_y_local <- sapply(local_depths, function(y){mean(local_depths_samp <= y)})



plot(points, local_depths, type="l")
points(points, predict(lcd_smooth), type="l", col="red")


depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:5000]),
                              beta = 0.5)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.05)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})

plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})


sum(depths_dens * r_y_dens)*0.02
sum(depths_dens * r_y_local)*0.02
sum(depths_dens * r_y_pvb)*0.02
sum(depths_dens * r_y_md)*0.02

sum(-1*derivs * r_y_local)*0.02 # 0.318
sum(-1*derivs * r_y_dens)*0.02 # 0.318
sum(-1*derivs * r_y_pvb)*0.02 # 0.318
sum(-1*derivs * r_y_md)*0.02 # 0.318

sum(-1*derivs_3 * r_y_local)*0.02 # 0
sum(-1*derivs_3 * r_y_dens)*0.02 # 0
sum(-1*derivs_3 * r_y_pvb)*0.02 # 0
sum(-1*derivs_3 * r_y_md)*0.02 # 0



p1_1 <- data.frame(x = points, y = depths_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Density", title = "Univariate Gaussian") +
  theme(text = element_text(size = 15))

p1_2 <- data.frame(x = points, y = local_depths) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.1, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p1_3 <- data.frame(x = points, y = md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p1_4 <- data.frame(x = points, y = depths_pvb) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.12, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "PVB local depth") +
  theme(text = element_text(size = 15))

p1_5 <- data.frame(x = points, y = r_y_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Density as depth") +
  theme(text = element_text(size = 15))

p1_6 <- data.frame(x = points, y = r_y_local) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.05, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p1_7 <- data.frame(x = points, y = r_y_md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p1_8 <- data.frame(x = points, y = r_y_pvb) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.07, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "PVB local depth") +
  theme(text = element_text(size = 15))

p1_deriv <- data.frame(x = points, y = derivs) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Derivative", title = "Derivative wrt scale parameter, Scenario 1") +
  theme(text = element_text(size = 15))

#pdf("example_1_1.pdf", width=18, height=6)
#grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4)
#dev.off()







### Scenario 2



mu <- c(0, 5)
alpha <- c(0.5, 0.5)
sd <- c(1, 1)


points <- seq(-5, 10, 0.02)
derivs <- deriv_fun(points)
derivs_2 <- deriv_fun_2(points)
derivs_3 <- deriv_fun_3(points)


nsamp = 20000
samp1 <- data_fun(nsamp)
samp2 <- data_fun(nsamp)
samp_points <- data_fun(nsamp)

depths_dens <- density_fun(points)
depths_samp_dens <- density_fun(samp_points)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})


local_depths <- rep(NA, length(points))
for(i in 1:length(points)){
  local_depths[i] <- local_depth_fun(points[i])
  print(i)
}
lcd_smooth <- loess(local_depths ~ points, span = 0.05)
local_depths_samp <- predict(lcd_smooth, samp_points)
r_y_local <- sapply(local_depths, function(y){mean(local_depths_samp <= y)})



plot(points, local_depths, type="l")
points(points, predict(lcd_smooth), type="l", col="red")


depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:5000]),
                              beta = 0.5)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.01)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})

plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})


sum(depths_dens * r_y_dens)*0.02
sum(depths_dens * r_y_local)*0.02
sum(depths_dens * r_y_pvb)*0.02
sum(depths_dens * r_y_md)*0.02

sum(-1*derivs * r_y_local)*0.02 # 0.304
sum(-1*derivs * r_y_dens)*0.02 # 0.316
sum(-1*derivs * r_y_pvb)*0.02 # 0.248
sum(-1*derivs * r_y_md)*0.02 # 0

sum(-1*derivs_2 * r_y_local)*0.02 # 0
sum(-1*derivs_2 * r_y_dens)*0.02 # 0
sum(-1*derivs_2 * r_y_pvb)*0.02 # 0
sum(-1*derivs_2 * r_y_md)*0.02 # 0

sum(-1*derivs_3 * r_y_local)*0.02 # 0
sum(-1*derivs_3 * r_y_dens)*0.02 # 0
sum(-1*derivs_3 * r_y_pvb)*0.02 # 0
sum(-1*derivs_3 * r_y_md)*0.02 # 0



p2_1 <- data.frame(x = points, y = depths_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Density", title = "Gaussian mixture") +
  theme(text = element_text(size = 15))

p2_2 <- data.frame(x = points, y = local_depths) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.1, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p2_3 <- data.frame(x = points, y = md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p2_4 <- data.frame(x = points, y = predict(pvb_smooth)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  #geom_smooth(se=F, span = 0.01, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "PVB local depth") +
  theme(text = element_text(size = 15))

p2_5 <- data.frame(x = points, y = r_y_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Density as depth") +
  theme(text = element_text(size = 15))

p2_6 <- data.frame(x = points, y = r_y_local) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.05, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p2_7 <- data.frame(x = points, y = r_y_md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p2_8 <- data.frame(x = points, y = r_y_pvb) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_line(lwd=1) +
  geom_smooth(se=F, span = 0.05, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "PVB local depth") +
  theme(text = element_text(size = 15))

p2_deriv <- data.frame(x = points, y = derivs) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Derivative", title = "Derivative wrt scale parameter, Scenario 2") +
  theme(text = element_text(size = 15))



depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:10000]),
                              beta = 0.1)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.05)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})
plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

sum(-1*derivs * r_y_pvb)*0.02 # 0.296








# Scenario 3


mu <- c(0, 5)
alpha <- c(0.5, 0.5)
sd <- c(1, 0.25)


points <- seq(-5, 10, 0.02)
derivs <- deriv_fun(points)
derivs_2 <- deriv_fun_2(points)
derivs_3 <- deriv_fun_3(points)


nsamp = 20000
samp1 <- data_fun(nsamp)
samp2 <- data_fun(nsamp)
samp_points <- data_fun(nsamp)

depths_dens <- density_fun(points)
depths_samp_dens <- density_fun(samp_points)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})


local_depths <- rep(NA, length(points))
for(i in 1:length(points)){
  local_depths[i] <- local_depth_fun(points[i])
  print(i)
}
lcd_smooth <- loess(local_depths ~ points, span = 0.05)
local_depths_samp <- predict(lcd_smooth, samp_points)
r_y_local <- sapply(local_depths, function(y){mean(local_depths_samp <= y)})



plot(points, local_depths, type="l")
points(points, predict(lcd_smooth), type="l", col="red")


depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:5000]),
                              beta = 0.5)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.01)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})

plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})


sum(depths_dens * r_y_dens)*0.02
sum(depths_dens * r_y_local)*0.02
sum(depths_dens * r_y_pvb)*0.02
sum(depths_dens * r_y_md)*0.02

sum(-1*derivs * r_y_local)*0.02 # 0.323
sum(-1*derivs * r_y_dens)*0.02 # 0.232
sum(-1*derivs * r_y_pvb)*0.02 # 0.309
sum(-1*derivs * r_y_md)*0.02 # 0


sum(-1*derivs_2 * r_y_local)*0.02 # 0.06
sum(-1*derivs_2 * r_y_dens)*0.02 # 0.436
sum(-1*derivs_2 * r_y_pvb)*0.02 # 0
sum(-1*derivs_2 * r_y_md)*0.02 # 0

sum(-1*derivs_3 * r_y_local)*0.02 # 0
sum(-1*derivs_3 * r_y_dens)*0.02 # 0
sum(-1*derivs_3 * r_y_pvb)*0.02 # 0
sum(-1*derivs_3 * r_y_md)*0.02 # 0.21



p3_1 <- data.frame(x = points, y = depths_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Density", title = "Gaussian mixture") +
  theme(text = element_text(size = 15))

p3_2 <- data.frame(x = points, y = local_depths) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  #geom_smooth(se=F, span = 0.1, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p3_3 <- data.frame(x = points, y = md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p3_4 <- data.frame(x = points, y = depths_pvb) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  #geom_smooth(se=F, span = 0.12, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "D(y,F)", title = "PVB local depth") +
  theme(text = element_text(size = 15))

p3_5 <- data.frame(x = points, y = r_y_dens) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Density as depth") +
  theme(text = element_text(size = 15))

p3_6 <- data.frame(x = points, y = r_y_local) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  #geom_smooth(se=F, span = 0.05, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Local community depth") +
  theme(text = element_text(size = 15))

p3_7 <- data.frame(x = points, y = r_y_md) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "Mahalanobis depth") +
  theme(text = element_text(size = 15))

p3_8 <- data.frame(x = points, y = r_y_pvb) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  #geom_smooth(se=F, span = 0.07, color="black", lwd=1) +
  theme_classic() +
  labs(x = "y", y = "R(y,F)", title = "PVB local depth") +
  theme(text = element_text(size = 15))

p3_deriv <- data.frame(x = points, y = derivs) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(lwd=1) +
  theme_classic() +
  labs(x = "y", y = "Derivative", title = "Derivative wrt scale parameter") +
  theme(text = element_text(size = 15))


grid.arrange(p3_1, p3_2, p3_3, p3_4, p3_5, p3_6, p3_7, p3_8, ncol=4)

library(patchwork)
pdf(file = "example_1_expanded.pdf", width=22, height = 28)
((p1_1 + labs(tag = "(a)")) | p1_2 | p1_3 | p1_4)/
  (p1_5 | p1_6 | p1_7 | p1_8)/
  ((p2_1 + labs(tag = "(b)")) | p2_2 | p2_3 | p2_4)/
  (p2_5 | p2_6 | p2_7 | p2_8)/
  ((p3_1 + labs(tag = "(c)")) | p3_2 | p3_3 | p3_4)/
  (p3_5 | p3_6 | p3_7 | p3_8) /
  ((p1_deriv + labs(tag = "(d)")) | p2_deriv | (p3_deriv + 
                                                  labs(title = "Derivative wrt scale parameter, Scenario 3")))
dev.off()





# Scenario 4


mu <- c(0, 5)
alpha <- c(0.25, 0.75)
sd <- c(1, 1)


points <- seq(-5, 10, 0.02)
derivs <- deriv_fun(points)
derivs_2 <- deriv_fun_2(points)
derivs_3 <- deriv_fun_3(points)


nsamp = 20000
samp1 <- data_fun(nsamp)
samp2 <- data_fun(nsamp)
samp_points <- data_fun(nsamp)

depths_dens <- density_fun(points)
depths_samp_dens <- density_fun(samp_points)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})


local_depths <- rep(NA, length(points))
for(i in 1:length(points)){
  local_depths[i] <- local_depth_fun(points[i])
  print(i)
}
lcd_smooth <- loess(local_depths ~ points, span = 0.05)
local_depths_samp <- predict(lcd_smooth, samp_points)
r_y_local <- sapply(local_depths, function(y){mean(local_depths_samp <= y)})



plot(points, local_depths, type="l")
points(points, predict(lcd_smooth), type="l", col="red")


depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:5000]),
                              beta = 0.5)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.01)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})

plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})


sum(depths_dens * r_y_dens)*0.02
sum(depths_dens * r_y_local)*0.02
sum(depths_dens * r_y_pvb)*0.02
sum(depths_dens * r_y_md)*0.02

sum(-1*derivs * r_y_local)*0.02 # 0.208
sum(-1*derivs * r_y_dens)*0.02 # 0.266
sum(-1*derivs * r_y_pvb)*0.02 # 0.278
sum(-1*derivs * r_y_md)*0.02 # 0.038

sum(-1*derivs_2 * r_y_local)*0.02 # 0.480
sum(-1*derivs_2 * r_y_dens)*0.02 # 0.409
sum(-1*derivs_2 * r_y_pvb)*0.02 # 0.190
sum(-1*derivs_2 * r_y_md)*0.02 # 0.461


sum(-1*derivs_3 * r_y_local)*0.02 # 0
sum(-1*derivs_3 * r_y_dens)*0.02 # 0
sum(-1*derivs_3 * r_y_pvb)*0.02 # 0
sum(-1*derivs_3 * r_y_md)*0.02 # 0.128





# Scenario 5


mu <- c(0, 5)
alpha <- c(0.25, 0.75)
sd <- c(1, 0.25)


points <- seq(-5, 10, 0.02)
derivs <- deriv_fun(points)
derivs_2 <- deriv_fun_2(points)
derivs_3 <- deriv_fun_3(points)


nsamp = 20000
samp1 <- data_fun(nsamp)
samp2 <- data_fun(nsamp)
samp_points <- data_fun(nsamp)

depths_dens <- density_fun(points)
depths_samp_dens <- density_fun(samp_points)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})


local_depths <- rep(NA, length(points))
for(i in 1:length(points)){
  local_depths[i] <- local_depth_fun(points[i])
  print(i)
}
lcd_smooth <- loess(local_depths ~ points, span = 0.05)
local_depths_samp <- predict(lcd_smooth, samp_points)
r_y_local <- sapply(local_depths, function(y){mean(local_depths_samp <= y)})



plot(points, local_depths, type="l")
points(points, predict(lcd_smooth), type="l", col="red")


depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:5000]),
                              beta = 0.5)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.01)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})

plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})


sum(depths_dens * r_y_dens)*0.02
sum(depths_dens * r_y_local)*0.02
sum(depths_dens * r_y_pvb)*0.02
sum(depths_dens * r_y_md)*0.02

sum(-1*derivs * r_y_local)*0.02 # 0.228
sum(-1*derivs * r_y_dens)*0.02 # 0.221
sum(-1*derivs * r_y_pvb)*0.02 # 0.186
sum(-1*derivs * r_y_md)*0.02 # 0

sum(-1*derivs_2 * r_y_local)*0.02 # 0.49
sum(-1*derivs_2 * r_y_dens)*0.02 # 0.48
sum(-1*derivs_2 * r_y_pvb)*0.02 # 0.497
sum(-1*derivs_2 * r_y_md)*0.02 # 0.49


sum(-1*derivs_3 * r_y_local)*0.02 # 0
sum(-1*derivs_3 * r_y_dens)*0.02 # 0
sum(-1*derivs_3 * r_y_pvb)*0.02 # 0
sum(-1*derivs_3 * r_y_md)*0.02 # 0.614





# Scenario 6


mu <- c(0, 5)
alpha <- c(0.25, 0.75)
sd <- c(0.25, 1)


points <- seq(-5, 10, 0.02)
derivs <- deriv_fun(points)
derivs_2 <- deriv_fun_2(points)
derivs_3 <- deriv_fun_3(points)


nsamp = 20000
samp1 <- data_fun(nsamp)
samp2 <- data_fun(nsamp)
samp_points <- data_fun(nsamp)

depths_dens <- density_fun(points)
depths_samp_dens <- density_fun(samp_points)
r_y_dens <- sapply(depths_dens, function(y){mean(depths_samp_dens <= y)})


local_depths <- rep(NA, length(points))
for(i in 1:length(points)){
  local_depths[i] <- local_depth_fun(points[i])
  print(i)
}
lcd_smooth <- loess(local_depths ~ points, span = 0.05)
local_depths_samp <- predict(lcd_smooth, samp_points)
r_y_local <- sapply(local_depths, function(y){mean(local_depths_samp <= y)})



plot(points, local_depths, type="l")
points(points, predict(lcd_smooth), type="l", col="red")


depths_pvb <- rep(NA, length(points))
for(i in 1:length(points)){
  depths_pvb[i] <- depthLocal(as.matrix(points[i]), as.matrix(samp1[1:5000]),
                              beta = 0.5)
  print(i)
}
pvb_smooth <- loess(depths_pvb ~ points, span = 0.01)
depths_pvb_samp <- predict(pvb_smooth, samp_points)
r_y_pvb <- sapply(depths_pvb, function(y){mean(depths_pvb_samp <= y)})

plot(points, depths_pvb, type="l")
points(points, predict(pvb_smooth), type="l", col="red")

md <- mdepth.MhD(points, samp1)$dep
md_samp <- mdepth.MhD(samp_points, samp1)$dep
r_y_md <- sapply(md, function(y){mean(md_samp <= y)})


sum(depths_dens * r_y_dens)*0.02
sum(depths_dens * r_y_local)*0.02
sum(depths_dens * r_y_pvb)*0.02
sum(depths_dens * r_y_md)*0.02

sum(-1*derivs * r_y_local)*0.02 # 0.224
sum(-1*derivs * r_y_dens)*0.02 # 0.308
sum(-1*derivs * r_y_pvb)*0.02 # 0.253
sum(-1*derivs * r_y_md)*0.02 # 0

sum(-1*derivs_2 * r_y_local)*0.02 # 0.471
sum(-1*derivs_2 * r_y_dens)*0.02 # 0.22
sum(-1*derivs_2 * r_y_pvb)*0.02 # 0.149
sum(-1*derivs_2 * r_y_md)*0.02 # 0.495


sum(-1*derivs_3 * r_y_local)*0.02 # 0
sum(-1*derivs_3 * r_y_dens)*0.02 # 0
sum(-1*derivs_3 * r_y_pvb)*0.02 # 0
sum(-1*derivs_3 * r_y_md)*0.02 # 0.075






