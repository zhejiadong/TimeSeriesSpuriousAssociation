library(tidyr) #Allows for us to manipulate the data structure
library(data.table) #Allows for us to manipulate the data structure
library(ggplot2) #this makes better looking plots in R
library(dplyr)
library(haven)
library(astsa)
library(TTR)
library(forecast)
library(broom.mixed)
library(cowplot)
library(knitr)
library(TSA)
library(energy)
library(plotrix)


# Code for AR (2) model

## sample
set.seed(123456)
popn <- 1000 # population number
ar.par <- c(0.9, -0.9) # ar parameter
p <- length(ar.par) # parameter length
eps <- rnorm( popn + p, 0, 1) # epsilon

ts <- rnorm(popn+p, 0, 1) # generate time series

for(i in (p+1) : (popn+p)){
  ts[i] <-   sum(ar.par * ts[(i-1):(i-p)]) + eps[i] 
}

ts.plot(ts[(p+1):(popn + p)]) # plot time series


# AR(2)
set.seed(123456)
popn <- 1000
n.sim <- 500
ar.par <- c(0.7,0.35)
p <- length(ar.par)


ar.process <- function(ar.par, popn, mu = 0, sd = 1){
  p <- length(ar.par)
  eps.x <- rnorm((popn + p), mean = mu, sd = sd)
  eps.y <- rnorm((popn + p), mean = mu, sd = sd)
  x <- rnorm(popn+p, 0, 1)
  y <- rnorm(popn+p, 0, 1)
  for(i in (1+p):(popn+p)){
    x[i] <-  sum(ar.par * x[(i-1):(i-p)]) + eps.x[i]  
    y[i] <-  sum(ar.par * y[(i-1):(i-p)])  + eps.y[i]
  }
  x <- x[(1+p) : (popn + p)]
  y <- y[(1+p) : (popn + p)]
  
  return(list("x" = x, "y" = y, "ar.par" = ar.par,
              "mu" = mu, "sd" = sd))
}

x.ts <- matrix(0, nrow = n.sim, ncol = popn) # per row is one replication
y.ts <- matrix(0, nrow = n.sim, ncol = popn)

for(i in 1: n.sim){
  per.process <- ar.process(ar.par, popn)
  x.ts[i,] <- per.process$x
  y.ts[i,] <- per.process$y
}


cor.coef <- c()
dis.core <- c()
dis.p <- c()
lm.p <- c()
lm.beta <- c()
z.beta <- c()


for(i in 1:n.sim){
  cor.coef[i] <- cor(x.ts[i,], y.ts[i,])
  dis.core[i] <- dcor(x.ts[i,], y.ts[i,])
  dd <- dcorT.test(x.ts[i,], y.ts[i,])
  dis.p[i] <-dd$p.value
  lt <- lm(y.ts[i,] ~ x.ts[i,])
  ls <- summary(lt)
  lm.beta[i] <- ls$coefficients[2,"Estimate"]
  lm.p[i] <- ls$coefficients[2,4] ## p value for linear regression test, p value.
  z.beta[i] <- lm.beta[i]/ (ls$coefficients[2,"Std. Error"] * sd(x.ts[i,]))
}




# AR(0.7, 0.3)
pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.7,0.3).corr.pdf",
    width = 6,
    height = 6)
hist(cor.coef, main = "Correlation coefficients for AR(0.7, 0.3)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.7,0.3).dist.pdf",
    width = 6,
    height = 6)
hist(dis.core, main = "Distance correlations for AR(0.7, 0.3)",
     xlab = "Distance correlation", col = "lightblue", prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.7,0.3).beta.pdf",
    width = 6,
    height = 6)
hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(0.7, 0.3)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()


data.frame( "linear test" = mean(lm.p < 0.05),
            "dis test" = mean(dis.p < 0.05)) 


pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(-0.7,-0.3).corr.pdf",
    width = 6,
    height = 6)
hist(cor.coef, main = "Correlation coefficients for AR(-0.7, -0.3)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(-0.7,-0.3).dist.pdf",
    width = 6,
    height = 6)
hist(dis.core, main = "Distance correlations for AR(-0.7, -0.3)",
     xlab = "Distance correlation", col = "lightblue", prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(-0.7,-0.3).beta.pdf",
    width = 6,
    height = 6)
hist(z.beta, main = expression(bold(paste("Standardized ", beta,
                                          " coefficients for AR(-0.7, -0.3)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T, cex.lab = 1.3, 
     cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()




# AR(0.6, 0.4)

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.6,0.4).corr.pdf",
    width = 6,
    height = 6)
hist(cor.coef, main = "Correlation coefficients for AR(0.6, 0.4)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.6,0.4).dist.pdf",
    width = 6,
    height = 6)
hist(dis.core, main = "Distance correlations for AR(0.6, 0.4)",
     xlab = "Distance correlation", col = "lightblue", prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.6,0.4).beta.pdf",
    width = 6,
    height = 6)
hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(0.6, 0.4)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

# AR(-0.6, -0.4)

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(-0.6,-0.4).beta.pdf",
    width = 6,
    height = 6)
hist(cor.coef, main = "Correlation coefficients for AR(-0.6, -0.4)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.6,0.4).dist.pdf",
    width = 6,
    height = 6)
hist(dis.core, main = "Distance correlations for AR(-0.6, -0.4)",
     xlab = "Distance correlation", col = "lightblue", prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

pdf(file = "~/Dropbox (Brown)/2023 Spring RAship/Spring23_RA/Pictures/AR(0.6,0.4).beta.pdf",
    width = 6,
    height = 6)
hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(-0.6, -0.4)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T,
     cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.3, cex.sub = 1.5)
dev.off()

# AR(0.6, 0.3)
hist(cor.coef, main = "Correlation coefficients for AR(0.6, 0.3)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T)

hist(dis.core, main = "Distance correlations for AR(0.6, 0.3)",
     xlab = "Distance correlation", col = "lightblue", prob=T)

hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(0.6, 0.3)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T)

# AR(-0.6, -0.3)
hist(cor.coef, main = "Correlation coefficients for AR(-0.6, -0.3)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T)

hist(dis.core, main = "Distance correlations for AR(-0.6, -0.3)",
     xlab = "Distance correlation", col = "lightblue", prob=T)

hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(-0.6, -0.3)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T)


# AR (0.5, 0.3)

ts.plot(x.ts[1,])
hist(cor.coef, main = "Correlation coefficients for AR(0.5, 0.3)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T)

hist(dis.core, main = "Distance correlations for AR(0.5, 0.3)",
     xlab = "Distance correlation", col = "lightblue", prob=T)

hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(0.5, 0.3)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T)


# AR(-0.5, -0.3)
ts.plot(x.ts[1,])
hist(cor.coef, main = "Correlation coefficients for AR(-0.5, -0.3)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T)

hist(dis.core, main = "Distance correlations for AR(-0.5, -0.3)",
     xlab = "Distance correlation", col = "lightblue", prob=T)

hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(-0.5, -0.3)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T)


# AR (0.9, -0.9)
hist(cor.coef, main = "Correlation coefficients for AR(0.9, -0.9)",
     xlab = "Correlation coefficient", col = "lightblue",prob=T)

hist(dis.core, main = "Distance correlations for AR(0.9, -0.9)",
     xlab = "Distance correlation", col = "lightblue", prob=T)

hist(z.beta, main = expression(bold(paste("Standardized ", beta, " coefficients for AR(0.9, -0.9)"))),
     xlab = expression(paste("Standardized ", beta, " coefficients")) ,
     col = "lightblue",prob=T)

