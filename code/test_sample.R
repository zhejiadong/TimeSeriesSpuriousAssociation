library("TTR")
library("forecast")
library("forecast")
kings <- scan("http://robjhyndman.com/tsdldata/misc/kings.dat",skip=3)
kingstimeseries <- ts(kings)

kingstimeseries
births <- scan("http://robjhyndman.com/tsdldata/data/nybirths.dat")
birthstimeseries <- ts(births, frequency=12, start=c(1946,1))

souvenir <- scan("http://robjhyndman.com/tsdldata/data/fancy.dat")
souvenirtimeseries <- ts(souvenir, frequency=12, start=c(1987,1))
plot.ts(kingstimeseries)
plot.ts(birthstimeseries)
plot.ts(souvenirtimeseries)
logsouvenirtimeseries <- log(souvenirtimeseries)
plot.ts(logsouvenirtimeseries)

kingstimeseriesSMA5 <- SMA(kingstimeseries,n=5)
kingstimeseriesSMA3 <- SMA(kingstimeseries,n=3)
plot.ts(kingstimeseriesSMA5)
plot.ts(kingstimeseriesSMA3)

birthstimeseriescomponents <- decompose(birthstimeseries)

birthstimeseriescomponents$seasonal
plot(birthstimeseriescomponents)

birthstimeseriesseasonallyadjusted <- 
  birthstimeseries - birthstimeseriescomponents$seasonal
plot(birthstimeseriesseasonallyadjusted)


rain <- scan("http://robjhyndman.com/tsdldata/hurst/precip1.dat",skip=1)
rainseries <- ts(rain,start=c(1813))
plot.ts(rainseries)

rainseriesforecasts <- HoltWinters(rainseries, beta=FALSE, gamma=FALSE)
plot(rainseriesforecasts)



HoltWinters(rainseries, beta=FALSE, gamma=FALSE, l.start=23.56)

rainseriesforecasts2 <- forecast.HoltWinters(rainseriesforecasts, h=8)




# ARIMA
skirts <- scan("http://robjhyndman.com/tsdldata/roberts/skirts.dat",skip=5)
skirtsseries <- ts(skirts,start=c(1866))
plot.ts(skirtsseries)
skirtsseriesdiff1 <- diff(skirtsseries, differences=1)
plot.ts(skirtsseriesdiff1)

skirtsseriesdiff2 <- diff(skirtsseries, differences=2)
plot.ts(skirtsseriesdiff2)
## find d
kingtimeseriesdiff1 <- diff(kingstimeseries, differences=1)
plot.ts(kingtimeseriesdiff1)
## find p,q
acf(kingtimeseriesdiff1, lag.max=20)
acf(kingtimeseriesdiff1, lag.max=20, plot=FALSE)
pacf(kingtimeseriesdiff1, lag.max=20)


pacf(kingtimeseriesdiff1, lag.max=20, plot=FALSE)


auto.arima(kings)


volcanodust <- scan("http://robjhyndman.com/tsdldata/annual/dvi.dat", skip=1)

kingstimeseriesarima <- arima(kingstimeseries, order=c(0,1,1))
kingstimeseriesarima

kingstimeseriesforecasts <- forecast(kingstimeseriesarima, h = 5)


