# Chapter 11, part 1 R code


## Sample CCF of the simulated data, 
## where X leads Y by 2 time units
## and the errors are white noise

# Exhibit 11.11
#win.graph(width=4.875, height=2.5,pointsize=8)
set.seed(12345)
X=rnorm(105)
Y=zlag(X,2)+.5*rnorm(105)
X=ts(X[-(1:5)],start=1,freq=1)
Y=ts(Y[-(1:5)],start=1,freq=1)
ccf(X,Y,ylab='CCF')


# Exhibit 11.12
phi=seq(0,.95,.15)
rejection=2*(1-pnorm(1.96*sqrt((1-phi^2)/(1+phi^2))))
M=signif(rbind(phi,rejection),2)
rownames(M)=c("phi", "Error Rate")
M


# Exhibit 11.13
set.seed(23457)
correlation.v=NULL
B=1000
n=500
for (i in 1:B) {x=cumsum(arima.sim(model=list(ma=.8),n=n))
y=cumsum(arima.sim(model=list(ma=.8),n=n))
correlation.v=c(correlation.v,ccf(x,y,lag.max=1,plot=F)$acf[2])
}
hist(correlation.v,prob=T,xlab=expression(r[0](X,Y)))


# Milk/electricity example:  Spurious correlation


# Exhibit 11.14
data(milk)
data(electricity)
milk.electricity=ts.intersect(milk,log(electricity))
# The ts.intersect function merges serveral time series into a panel of time
# series over the time frame where each series has data.

plot(milk.electricity,yax.flip=T) 
# the option yax.flip=T flips the y-axis for the series alternately so as
# to make the labeling clearer.

# Exhibit 11.15
ccf(as.numeric(milk.electricity[,1]),as.numeric(milk.electricity[,2]),
    main='milk & electricity',ylab='CCF')
# The as.numeric function strips the time series attribute from the time series.
# This is done to nullify the default way the ccf function plots the cross-correlations.
# You may want to repeat the command without the as.numeric function to see 
# the default labels of the lags according to the period of the data.
# ccf((milk.electricity[,1]),(milk.electricity[,2]), main='milk & electricity',ylab='CCF')

## prewhitening to deal with the spurious correlation:

# Exhibit 11.16
me.dif=ts.intersect(diff(diff(milk,12)),diff(diff(log(electricity),12)))
prewhiten(as.numeric(me.dif[,1]),as.numeric(me.dif[,2]),
          ylab='CCF' )


## Bluebird Potato Chip example

# Exhibit 11.17
data(bluebird)
plot(bluebird,yax.flip=T)

sales=bluebird[,1]
price=bluebird[,2]

# Exhibit 11.18
prewhiten(y=diff(bluebird)[,1],x=diff(bluebird)[,2],main="Price & Sales",ylab='CCF')
# As the time series are of unit period, there is no need to apply the as.numeric 
# function.

# Note the 'sales' vector has ALREADY been log-transformed.

# Try: hist(sales)
# and: qqnorm(sales)
# to see that the transformation has achieved approximate symmetry/normality.

# Exhibit 11.19
sales=bluebird[,1]
price=bluebird[,2]
chip.m1=lm(sales~price,data=bluebird)
summary(chip.m1)

# Exhibit 11.20
acf(residuals(chip.m1),ci.type='ma')

# Exhibit 11.21
pacf(residuals(chip.m1))

# Exhibit 11.22
eacf(residuals(chip.m1))

# Exhibit 11.23
chip.m2=arima(sales,order=c(1,0,4),xreg=data.frame(price))
chip.m2
# MA(1)& MA(3) estimates are not significant, at 5% level,
# so they are constrained to be zero in the model fit below.

chip.m3=arima(sales,order=c(1,0,4),xreg=data.frame(price),fixed=c(NA,0,NA,0,NA,NA,NA))
# The MA(1) & MA(3) estimates are the second and fourth coefficients listed
# in the model fit chip.m2. They are set to be zero using the fixed option.
# NAs in the fixed option correspond to free parameters. 
chip.m3

# Now, the AR(1) coefficient estimate also seems not significant, so it is
# removed in the next fitted model.
chip.m4=arima(sales,order=c(0,0,4),xreg=data.frame(price),fixed=c(0,NA,0,NA,NA,NA))
chip.m4

# model diagnostics can be done by running the tsdiag command.
tsdiag(chip.m4)

# Forecasting the next 3 months of logged sales, 
# if the next 3 months' prices are 1.75, 1.70, and 1.68:

plot(n1=95,chip.m4,n.ahead=3,newxreg=data.frame(c(1.75,1.70, 1.68)))  # the plot of the forecasts
# n1=95 tells R to only plot the orginal time series values starting from time 95

# the actual forecasts and the standard errors:
predict(chip.m4,n.ahead=3,newxreg=data.frame(c(1.75,1.70, 1.68))) 


## A Shumway and Stoffer example:

library(astsa)
data(cmort)

trend = time(cmort); temp = tempr - mean(tempr); temp2 = temp^2
fit = lm(cmort~trend + temp + temp2 + part, na.action=NULL)  #OLS fit
acf2(resid(fit), 52) # implies AR2 process for the noise
#sarima(cmort, 2,0,0, xreg=cbind(trend, temp, temp2, part) )  # WLS fit using 'sarima'

fit.wls=arima(cmort,order=c(2,0,0),xreg=cbind(trend, temp, temp2, part) )  # WLS fit using 'arima'
fit.wls

# model diagnostics can be done by running the tsdiag command.

plot(resid(fit.wls))
tsdiag(fit.wls)
qqnorm(resid(fit.wls))

# Do the residuals from the WLS fit look like white noise?


## An example with a Lagged Covariate and SARIMA errors:

lag2.plot(soi, rec, 8) # to see the appropriate choice of lag

# The strongest linear association seems to be at lag 6.

# But notice the linear trend is only there when SOI > 0, 
# so we use a dummy (indicator) variable for that in our model:

dummy = ifelse(soi<0, 0, 1)
fish = ts.intersect(rec, soiL6=lag(soi,-6), dL6=lag(dummy,-6), dframe=TRUE)
summary(fit <- lm(rec ~soiL6*dL6, data=fish, na.action=NULL))
attach(fish)
plot(resid(fit))
acf2(resid(fit)) # indicates AR(2)
(fit1 = arima(rec,order=c(2,0,0), xreg = cbind(soiL6, dL6, I(soiL6*dL6) )))
acf2(resid(fit1)) # indicates seasonal AR
intract = soiL6*dL6 # interaction term
# try seasonal AR order 1 ...
(fit2 = arima(rec, order=c(2,0,0), seasonal=list(order=c(1,0,0),period=12), xreg = cbind(soiL6, dL6, intract) ))
acf2(resid(fit2))

# ... and then seasonal AR order 2, which appears to be best
(fit3 = arima(rec, order=c(2,0,0), seasonal=list(order=c(2,0,0),period=12), xreg = cbind(soiL6, dL6, intract) ))
acf2(resid(fit3))

plot(rstandard(fit3))



# Chapter 11, part 2 R code

library(TSA)

# Exhibit 11.1

## Plot of the airmiles time series

#win.graph(width=4.875, height=2.5,pointsize=8)
data(airmiles)
plot(log(airmiles),ylab='Log(airmiles)',xlab='Year', )
points(y=log(airmiles),x=time(airmiles),pch=as.vector(season(airmiles)))

library(astsa)

acf2(window(log(airmiles),end=c(2001,8)), 48)
# Seasonality is clearly apparent, 
# and so is the nonstationarity from the initial time series plot.

# # After taking first differences and seasonal differences, we see a pattern indicating an MA(1) model

# Exhibit 11.5
acf(as.vector(diff(diff(window(log(airmiles),end=c(2001,8)),12))),lag.max=48)

# ACF values cutting off after lag 1 ...

pacf(as.vector(diff(diff(window(log(airmiles),end=c(2001,8)),12))),lag.max=48)

# PACF values tailing off?

# So try ARIMA(0,1,1)x(0,1,0)_12 model for the pre-intervention data:


air.pre=arima(window(log(airmiles),end=c(2001,8)),order=c(0,1,1),seasonal=list(order=c(0,1,0),period=12) )

tsdiag(air.pre)

# Note one major outlier in the residuals of the seasonal model.  We will discuss more later.



### Handling Outliers

# Exhibit 11.9 
set.seed(12345)
y=arima.sim(model=list(ar=.8,ma=.5),n.start=158,n=100)
y[10]
y[10]=10
y=ts(y,freq=1,start=1)
plot(y,type='o')
acf(y)
pacf(y)
eacf(y)
m1=arima(y,order=c(1,0,0))
m1
detectAO(m1)
detectAO(m1, robust=F) # non-robust version of the test has less power to detect outliers
detectIO(m1)
tsdiag(m1)  # note the major outlier

m2=arima(y,order=c(1,0,0),xreg=data.frame(AO=seq(y)==10))
m2
detectAO(m2)
detectIO(m2)
tsdiag(m2)

m3=arima(y,order=c(1,0,1),xreg=data.frame(AO=seq(y)==10))
detectAO(m3)
detectIO(m3)
tsdiag(m3)
m3

plot(y,type='b') 
arrows(30,6, 11,9.8, length=.1,angle=20)
text(34,6.2, "AO")


# Exhibit 11.10
library(TSA)
data(co2)
m1.co2=arima(co2,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
m1.co2
tsdiag(m1.co2) # note the one large standardized residual

detectAO(m1.co2)
detectIO(m1.co2)

time(co2)
time(co2)[57]

m4.co2=arimax(co2,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),
              io=c(57))
m4.co2

tsdiag(m4.co2)

# the signs of the MA coefficient estimates appear to be opposite to 
# those in the book, because R uses the plus convention in the MA parameterization.

## Back to the airmiles example:


air.m1.tentative=arimax(log(airmiles),order=c(0,1,1),seasonal=list(order=c(0,1,0),
                                                                   period=12),xtransf=data.frame(I911=1*(seq(airmiles)==69),
                                                                                                 I911=1*(seq(airmiles)==69)),
                        transfer=list(c(0,0),c(1,0)),method='ML')

# Transfer function components are incorporated by the xtransf and transfer
# arguments.
# Here, the transfer function consists of two parts omega0*P(t) and 
# omega1/(1-omega2*B)P(t) where the inputs of the two transfer
# functions are identical and equals the dummy variable that is 1 at September
# 2001 (the 69th data point) and zero otherwise.
# xtransf is a matrix whose columns are the input variables.
# transfer is a list consisting of the pair of (AR order, MA order) of each
# transfer function, which in this examples is (0,0) and (1,0).

air.m1.tentative

tsdiag(air.m1.tentative)

# Note the significant lag-12 value in the ACF plot for the residuals.
# We should add a seasonal MA(1) term to account for this:

detectAO(air.m1.tentative)

detectIO(air.m1.tentative)

# We also include the three additive outliers using
# dummy variables in the xreg argument:

# The detectAO function marks observation 25 as an AO, 
# but due to the seasonal differencing and first differencing, 
# we see that the observations 12 and 13 are really the ones that break from the pattern:

plot(log(airmiles),ylab='Log(airmiles)',xlab='Year', )
points(y=log(airmiles),x=time(airmiles),pch=as.vector(season(airmiles)))
points(1996.917, log(airmiles)[12],col='red',cex=2.5)
points(1997, log(airmiles)[13],col='red',cex=2.5)
points(2002.917, log(airmiles)[84],col='red',cex=2.5)

# Exhibit 11.6 
air.m1=arimax(log(airmiles),order=c(0,1,1),seasonal=list(order=c(0,1,1),
                                                         period=12),xtransf=data.frame(I911=1*(seq(airmiles)==69),
                                                                                       I911=1*(seq(airmiles)==69)),
              transfer=list(c(0,0),c(1,0)),xreg=data.frame(Dec96=1*(seq(airmiles)==12),
                                                           Jan97=1*(seq(airmiles)==13),Dec02=1*(seq(airmiles)==84)),method='ML')

# Additive outliers are incorporated as dummy variables in xreg.
# Transfer function components are incorporated by the xtransf and transfer
# arguments.
# Here, the transfer function consists of two parts omega0*P(t) and 
# omega1/(1-omega2*B)P(t) where the inputs of the two transfer
# functions are identical and equals the dummy variable that is 1 at September
# 2001 (the 69th data point) and zero otherwise.
# xtransf is a matrix whose columns are the input variables.
# transfer is a list consisting of the pair of (AR order, MA order) of each
# transfer function, which in this example is (0,0) and (1,0).

# Note in our model for m_t:
# The first term has order 0 for B in both numerator and denominator (B does not appear)
# The second term has order 0 for B in the numerator and 1 in the denominator.

# See pages 451-454 of the Cryer and Chan textbook for full details of specifying
# the transfer argument with an intervention model.

air.m1

# Exhibit 11.7
plot(log(airmiles),ylab="log(airmiles)")
# The line plot shows the observed data.
points(fitted(air.m1))
# The open circles are the fitted values under the final model.
# Looks like a fairly good fit.

## Plotting the intervention effect over time:

# Exhibit 11.8
Nine11p=1*(seq(airmiles)==69)
plot(ts(Nine11p*(-0.0949)+
          filter(Nine11p,filter=.8139,method='recursive',side=1)*(-0.2715),
        frequency=12,start=1996),type='h',ylab='9/11 Effects')
abline(h=0)


## A time series regression example with an outlier:


# Exhibit 11.24
data(boardings)
plot(boardings,yax.flip=T)
# The Denver dataset has three time series. Here, we only plot the first 
# two series.

log.boardings=boardings[,1]
log.price=boardings[,2]

library(astsa)

acf2(diff(log.price))

# ACF has damped sine wave, PACF cuts off after 2 lags
# AR(2) on first differences --> ARIMA(2,1,0) for logged prices.

# Exhibit 11.25

# Prewhitening based on ARIMA(2,1,0) for logged prices:
m1=arima(boardings[,2],order=c(2,1,0))
prewhiten(x=boardings[,2],y=boardings[,1],x.model=m1)

# Prewhitening based on generically fitting a long AR model:
prewhiten(x=diff(boardings)[,2],y=diff(boardings)[,1],ylab='CCF')


# Exhibit 11.26
log.boardings=boardings[,1]
log.price=boardings[,2]

# Initially just fitting an OLS model of log boardings on log price:

boardings.m.lm = lm(log.boardings~log.price)
acf2(residuals(boardings.m.lm))

# Early PACFs cut off after lag 2, but note significant PACF at lag 12 (and 13).
#  ---> seasonal $ARIMA(2,0,0) \times (1,0,0)_{12}$ model for the noise process.

boardings.m1.tent=arima(log.boardings,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=12),
                        xreg=data.frame(log.price))
boardings.m1.tent

# AR(2) coefficient not significant, so we remove it, and this improves the model:

boardings.m1=arima(log.boardings,order=c(1,0,0),seasonal=list(order=c(1,0,0),period=12),
                   xreg=data.frame(log.price))
boardings.m1

tsdiag(boardings.m1)
plot(rstandard(boardings.m1))  # One or more notable residuals

detectAO(boardings.m1) 
detectIO(boardings.m1)
# An AO is detected at time point 32, as well as an IO detected at time point 44
# Since the test statistic of the AO has larger magnitude, an AO is added to the
# model fitted below.

dummy.32 = c(rep(0,31),1,rep(0,36)) # takes value 1 only at time 32

boardings.m2.tent1=arima(log.boardings,order=c(1,0,0),seasonal=list(order=c(1,0,0),period=12),
                         xreg=data.frame(log.price,outlier=dummy.32))
boardings.m2.tent1

tsdiag(boardings.m2.tent1)

# Significant lag-3 ACF value of residuals, so we add an MA(3) component:

boardings.m2.tent2=arima(log.boardings,order=c(1,0,3),seasonal=list(order=c(1,0,0),period=12),
                         xreg=data.frame(log.price,outlier=dummy.32))
boardings.m2.tent2

tsdiag(boardings.m2.tent2)

# MA1 and MA2 coefficients not significant, so we can fix them at 0:

boardings.m2=arima(log.boardings,order=c(1,0,3),seasonal=list(order=c(1,0,0),period=12),
                   xreg=data.frame(log.price,outlier=dummy.32),
                   fixed=c(NA,0,0,rep(NA,5)))
boardings.m2
detectAO(boardings.m2)
detectIO(boardings.m2)
# No outliers are detected!
tsdiag(boardings.m2,tol=.15,gof.lag=24)
# Model diagnostics appear to be OK.

## Compare the model with the outlier unaccounted for:
boardings.m1

## to the model with the outlier incorporated into it:
boardings.m2

