# ODIN baseline and colocation analysis - June 2015
Gustavo Olivares  


```r
###############################################################################
# Analysis script for the manuscript "The Outdoor Dust Information Node (ODIN)
# - Development and performance assessment of a low cost ambient dust sensor" 
# submitted to Atmospheric Measurement Techniques. June 2015
###############################################################################
## Preamble ####
# Libraries
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
require(openair)
```

```
## Loading required package: openair
## Loading required package: lazyeval
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: maps
```

```r
require(reshape2)
```

```
## Loading required package: reshape2
```

```r
require(gridExtra)
```

```
## Loading required package: gridExtra
## Loading required package: grid
```

```r
# Seed
set.seed(1)
###############################################################################

## Baseline ####
raw_dust_data_baseline <- read.table("http://files.figshare.com/2099394/sharp_baseline_test.tsv",
                            header = T,
                            sep = "",
                            quote = "\"'",
                            dec = ".",
                            row.names=NULL)
# Scale Dust measurements to mV (digital scale is 0 - 1024 translates to 0 - 5000 mV)
for (i in (3:10)){
  raw_dust_data_baseline[,i]<-raw_dust_data_baseline[,i]*5000/1024
}
# force GMT as the time zone to avoid openair issues with daylight saving switches
# The actual time zone is 'NZST'
raw_dust_data_baseline$date=as.POSIXct(paste(raw_dust_data_baseline$Date,raw_dust_data_baseline$Time),tz='GMT')
raw_dust_data_baseline$Time<-NULL
raw_dust_data_baseline$Date<-NULL

# Averages for 10 minutes
dust_data_baseline<-timeAverage(mydata = raw_dust_data_baseline,avg.time = '10 min',statistic = 'mean')
dust_summary_baseline <- lapply( dust_data_baseline[,(2:9)] , function(x) rbind( mean = mean(x),
                                                               sd = sd(x) ,
                                                               median = median(x),
                                                               minimum = min(x),
                                                               maximum = max(x),
                                                               s.size = length(x)))
dust_summary_baseline <- data.frame( dust_summary_baseline )
dust_summary_baseline
```

```
##               Dust1       Dust2       Dust3       Dust4       Dust5
## mean      0.2429041   0.9865411  27.9243672   0.9608271   8.6837857
## sd        0.1525655   0.3008985   0.9117687   0.2750282   0.3684347
## median    0.2077793   0.9350066  27.9463098   0.9350066   8.7267287
## minimum   0.0000000   0.2077793  25.0374003   0.2077793   7.2722739
## maximum   0.7272274   1.7661237  30.4396609   1.9739029   9.6617354
## s.size  678.0000000 678.0000000 678.0000000 678.0000000 678.0000000
##               Dust6       Dust7       Dust8
## mean      4.4550723  10.9696956  10.6632976
## sd        0.3009904   0.3966578   0.4504248
## median    4.4672540  10.9084109  10.7006316
## minimum   3.5603841   9.7656250   9.3500665
## maximum   5.5061503  12.0511968  12.1053060
## s.size  678.0000000 678.0000000 678.0000000
```

```r
###############################################################################

## Temperature dependence ####
raw_dust_data_temperature <- read.table("http://files.figshare.com/2099396/sharp_baseline_temperature_Nov2014.tsv",
                            header = T,
                            sep = "",
                            quote = "\"'",
                            dec = ".",
                            row.names=NULL)
# Baseline estimates from previous tests
baseline<-as.numeric(dust_summary_baseline['mean',])
# Scale Dust measurements to mV (digital scale is 0 - 1024 translates to 0 - 5000 mV)
# Substract previous baseline data from the signal.
for (i in (3:10)){
  raw_dust_data_temperature[,i]<-raw_dust_data_temperature[,i]*5000/1024-baseline[i-2]
}
# force GMT as the time zone to avoid openair issues with daylight saving switches
# The actual time zone is 'NZST'
raw_dust_data_temperature$date=as.POSIXct(paste(raw_dust_data_temperature$Date,raw_dust_data_temperature$Time),tz='GMT')
raw_dust_data_temperature$Time<-NULL
raw_dust_data_temperature$Date<-NULL

# Averages for 10 minutes
dust_data_temperature<-timeAverage(mydata = raw_dust_data_temperature,avg.time = '10 min',statistic = 'mean')
long_scatter_temperature=reshape(dust_data_temperature[,c(2,3,4,5,6,7,8,9,10,11)],direction = 'long',varying = names(dust_data_temperature[2:9]),v.name='Dust',timevar = 'Sensor')
long_scatter_temperature$Sensor=as.character(long_scatter_temperature$Sensor)
png("./sensor_temperature_scatter.png",height = 620,width = 1280)
ggplot(long_scatter_temperature, aes(x=Temperature, y=Dust, color=Sensor)) +
  geom_point(shape=19,alpha=1/4)+
  theme(text = element_text(size = 30))+
  xlab('Temperature (C)')+
  ylab('Sensor response (mV)')
```

```
## Warning in loop_apply(n, do.ply): Removed 744 rows containing missing
## values (geom_point).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
summary(Temp.lm<-lm(data=long_scatter_temperature,Dust~Temperature))
```

```
## 
## Call:
## lm(formula = Dust ~ Temperature, data = long_scatter_temperature)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -6.4725 -2.9810 -0.4161  2.5320  9.4704 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  3.25219    0.18366   17.71   <2e-16 ***
## Temperature  0.27018    0.01145   23.61   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 3.462 on 2742 degrees of freedom
##   (744 observations deleted due to missingness)
## Multiple R-squared:  0.1689,	Adjusted R-squared:  0.1686 
## F-statistic: 557.3 on 1 and 2742 DF,  p-value: < 2.2e-16
```

```r
###############################################################################

## ODIN's performance evaluation ####
odin_01 <- read.table("http://files.figshare.com/2099373/odin01_24July_14August_2014.tsv",
                      header=TRUE,
                      quote="")
ecan_data<-read.table("http://files.figshare.com/2099539/StAlbans24July_14August_2014.tsv",
                      header = TRUE,
                      sep="\t",
                      quote = "\"'",
                      dec=".")
ecan_data$date<-as.POSIXct(ecan_data$date,format = "%m/%d/%Y %H:%M:%S")
odin_01$date=as.POSIXct(paste(odin_01$Date,odin_01$Time),tz='GMT')
odin_01$Time<-NULL
odin_01$Date<-NULL
odin_01$Battery<-5*odin_01$Battery/1024
# Merging the data
all_merged_FULL<-merge(odin_01,ecan_data,by = 'date',all = TRUE)
all_merged.10min<-timeAverage(all_merged_FULL,avg.time = '10 min')
# Time sync
# Check the time difference, correct the data and re-merge.
lag_test=ccf(all_merged.10min$Temperature,
             all_merged.10min$Temperature.1m,
             na.action=na.pass,
             lag.max=100,
             type='correlation',
             ylab='Correlation',
             main='Temperature correlation as function of clock lag')
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-1.png) 

```r
odin01_lag=lag_test$lag[which.max(lag_test$acf)]
odin_01$date=odin_01$date-odin01_lag*10*60
odin_01$ODIN = odin_01$Dust
odin_01$Dust<-NULL
all_merged_FULL<-merge(odin_01,ecan_data,by = 'date',all = TRUE)
all_merged.10min<-timeAverage(all_merged_FULL,avg.time = '10 min')

# Dust performance using ECan data for calibration
# Calibration expression:
#  $ODIN_{calibrated}=A*ODIN_{raw}+B*Temperature_{ODIN}+C*RH_{ODIN}+D$

# Full dataset 1 hour  
all_merged.1hr<-timeAverage(all_merged.10min,avg.time='1 hour')
all_merged.1hr$PMcoarse<-all_merged.1hr$PM10.FDMS - all_merged.1hr$PM2.5.FDMS
# PM$_{2.5}$ fdms
summary(odin1.lm.1hr.pm2.5.fdms<-
          lm(data=all_merged.1hr,PM2.5.FDMS~
               ODIN+Temperature+Humidity))
```

```
## 
## Call:
## lm(formula = PM2.5.FDMS ~ ODIN + Temperature + Humidity, data = all_merged.1hr)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -32.901  -7.623  -1.608   5.541  95.078 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -35.57353    7.55221  -4.710 3.21e-06 ***
## ODIN          0.37554    0.01526  24.605  < 2e-16 ***
## Temperature  -4.45436    0.17393 -25.610  < 2e-16 ***
## Humidity     -0.63004    0.07362  -8.558  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13.79 on 502 degrees of freedom
##   (22 observations deleted due to missingness)
## Multiple R-squared:  0.7002,	Adjusted R-squared:  0.6984 
## F-statistic: 390.8 on 3 and 502 DF,  p-value: < 2.2e-16
```

```r
format(coef(odin1.lm.1hr.pm2.5.fdms),digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     "-35.6"     "  0.4"     " -4.5"     " -0.6"
```

```r
format((confint(odin1.lm.1hr.pm2.5.fdms)[,2]-confint(odin1.lm.1hr.pm2.5.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     "14.84"     " 0.03"     " 0.34"     " 0.14"
```

```r
# PM$_{10}$ fdms
summary(odin1.lm.1hr.pm10.fdms<-
          lm(data=all_merged.1hr,PM10.FDMS~
               ODIN+Temperature+Humidity))
```

```
## 
## Call:
## lm(formula = PM10.FDMS ~ ODIN + Temperature + Humidity, data = all_merged.1hr)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -35.859 -10.585  -2.442   7.733  97.719 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -43.12180    9.69449  -4.448 1.07e-05 ***
## ODIN          0.42655    0.01959  21.771  < 2e-16 ***
## Temperature  -4.47796    0.22327 -20.056  < 2e-16 ***
## Humidity     -0.64501    0.09451  -6.825 2.54e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 17.7 on 502 degrees of freedom
##   (22 observations deleted due to missingness)
## Multiple R-squared:  0.6084,	Adjusted R-squared:  0.6061 
## F-statistic:   260 on 3 and 502 DF,  p-value: < 2.2e-16
```

```r
format(coef(odin1.lm.1hr.pm10.fdms),digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     "-43.1"     "  0.4"     " -4.5"     " -0.6"
```

```r
format((confint(odin1.lm.1hr.pm10.fdms)[,2]-confint(odin1.lm.1hr.pm10.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     "19.05"     " 0.04"     " 0.44"     " 0.19"
```

```r
# PM$_{coarse}$ fdms
summary(odin1.lm.1hr.pmcoarse.fdms<-
          lm(data=all_merged.1hr,PMcoarse~
               ODIN+Temperature+Humidity))
```

```
## 
## Call:
## lm(formula = PMcoarse ~ ODIN + Temperature + Humidity, data = all_merged.1hr)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -15.643  -4.900  -1.853   2.623  55.116 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.548270   4.487454  -1.682   0.0932 .  
## ODIN         0.051013   0.009069   5.625 3.09e-08 ***
## Temperature -0.023600   0.103349  -0.228   0.8195    
## Humidity    -0.014970   0.043746  -0.342   0.7323    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.193 on 502 degrees of freedom
##   (22 observations deleted due to missingness)
## Multiple R-squared:  0.07253,	Adjusted R-squared:  0.06699 
## F-statistic: 13.09 on 3 and 502 DF,  p-value: 3.064e-08
```

```r
format(coef(odin1.lm.1hr.pmcoarse.fdms),digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     "-7.55"     " 0.05"     "-0.02"     "-0.01"
```

```r
format((confint(odin1.lm.1hr.pmcoarse.fdms)[,2]-confint(odin1.lm.1hr.pmcoarse.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##      "8.82"      "0.02"      "0.20"      "0.09"
```

```r
# 24 hour dataset
merged.24hr <- timeAverage(all_merged.1hr,avg.time = '24 hour',statistic = 'mean')

# PM$_{2.5}$ fdms
summary(odin1.lm.1dy.pm2.5.fdms<-
          lm(data=merged.24hr,PM2.5.FDMS~
               ODIN+Temperature+Humidity))
```

```
## 
## Call:
## lm(formula = PM2.5.FDMS ~ ODIN + Temperature + Humidity, data = merged.24hr)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -13.2703  -4.0617  -0.2754   4.1605  17.2774 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -3.16933   31.38948  -0.101   0.9207    
## ODIN         0.31459    0.05786   5.437 3.64e-05 ***
## Temperature -4.60105    0.86697  -5.307 4.80e-05 ***
## Humidity    -0.79303    0.36717  -2.160   0.0445 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.715 on 18 degrees of freedom
## Multiple R-squared:  0.7081,	Adjusted R-squared:  0.6595 
## F-statistic: 14.56 on 3 and 18 DF,  p-value: 4.649e-05
```

```r
format(coef(odin1.lm.1dy.pm2.5.fdms),digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##      "-3.2"      " 0.3"      "-4.6"      "-0.8"
```

```r
format((confint(odin1.lm.1dy.pm2.5.fdms)[,2]-confint(odin1.lm.1dy.pm2.5.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##      "65.9"      " 0.1"      " 1.8"      " 0.8"
```

```r
# PM$_{10}$ fdms
summary(odin1.lm.1dy.pm10.fdms<-
          lm(data=merged.24hr,PM10.FDMS~
               ODIN+Temperature+Humidity))
```

```
## 
## Call:
## lm(formula = PM10.FDMS ~ ODIN + Temperature + Humidity, data = merged.24hr)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -17.2492  -5.5529   0.0465   4.7099  18.7966 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  2.58157   36.39163   0.071 0.944229    
## ODIN         0.35067    0.06708   5.228 5.69e-05 ***
## Temperature -4.76316    1.00513  -4.739 0.000164 ***
## Humidity    -0.92979    0.42568  -2.184 0.042420 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 10.1 on 18 degrees of freedom
## Multiple R-squared:  0.6679,	Adjusted R-squared:  0.6126 
## F-statistic: 12.07 on 3 and 18 DF,  p-value: 0.0001448
```

```r
format(coef(odin1.lm.1dy.pm10.fdms),digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##      " 2.6"      " 0.4"      "-4.8"      "-0.9"
```

```r
format((confint(odin1.lm.1dy.pm10.fdms)[,2]-confint(odin1.lm.1dy.pm10.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##      "76.5"      " 0.1"      " 2.1"      " 0.9"
```

```r
# PM$_{coarse}$ fdms
summary(odin1.lm.1dy.pmcoarse.fdms<-
          lm(data=merged.24hr,PMcoarse~
               ODIN+Temperature+Humidity))
```

```
## 
## Call:
## lm(formula = PMcoarse ~ ODIN + Temperature + Humidity, data = merged.24hr)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -4.4221 -2.4466 -0.6244  1.8651  9.0318 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  5.75090   14.21592   0.405    0.691
## ODIN         0.03608    0.02620   1.377    0.185
## Temperature -0.16212    0.39264  -0.413    0.685
## Humidity    -0.13675    0.16629  -0.822    0.422
## 
## Residual standard error: 3.947 on 18 degrees of freedom
## Multiple R-squared:  0.1776,	Adjusted R-squared:  0.04053 
## F-statistic: 1.296 on 3 and 18 DF,  p-value: 0.3064
```

```r
format(coef(odin1.lm.1dy.pmcoarse.fdms),digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     " 5.75"     " 0.04"     "-0.16"     "-0.14"
```

```r
format((confint(odin1.lm.1dy.pmcoarse.fdms)[,2]-confint(odin1.lm.1dy.pmcoarse.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN Temperature    Humidity 
##     "29.87"     " 0.06"     " 0.82"     " 0.35"
```

```r
### Calculate the corrected data
all_merged.1hr$ODIN.2.5f<-predict(odin1.lm.1hr.pm2.5.fdms,
                                  newdata = all_merged.1hr)
all_merged.1hr$ODIN.10f<-predict(odin1.lm.1hr.pm10.fdms,
                                 newdata = all_merged.1hr)
all_merged.1hr$ODIN.coarse<-predict(odin1.lm.1hr.pmcoarse.fdms,
                                    newdata = all_merged.1hr)

merged.24hr$ODIN.2.5f<-predict(odin1.lm.1dy.pm2.5.fdms,
                                  newdata = merged.24hr)
merged.24hr$ODIN.10f<-predict(odin1.lm.1dy.pm10.fdms,
                                 newdata = merged.24hr)
merged.24hr$ODIN.coarse<-predict(odin1.lm.1dy.pmcoarse.fdms,
                              newdata = merged.24hr)

png('./raw_odin_fdms.png',width = 1500,height = 750)
tseries_pm2.5<-ggplot(data = all_merged.1hr, aes(x=date,y=PM2.5.FDMS))+
  geom_line(colour = 'red',linetype = 2)+
  ylab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
tseries_ODIN<-ggplot(data = all_merged.1hr, aes(x=date,y=ODIN))+
  geom_line(colour = 'black',linetype = 1)+
  ylab('ODIN(mV)')+
  theme(text=element_text(size=30))
tseries_pm10<-ggplot(data = all_merged.1hr, aes(x=date,y=PM10.FDMS))+
  geom_line(colour = 'blue',linetype = 2)+
  ylab(bquote(PM[10]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
tseries_pmcoarse<-ggplot(data = all_merged.1hr, aes(x=date,y=PMcoarse))+
  geom_line(colour = 'blue',linetype = 2)+
  ylab(bquote(PM[coarse]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
grid.arrange(tseries_pm2.5,tseries_ODIN,tseries_pmcoarse,ncol = 1)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
png('./corrected_odin_fdms_PM2.5.png',width = 2400,height = 1200)
tseries_pm2.5_1hr<-ggplot(data = all_merged.1hr, aes(x=date))+
  geom_line(aes(y=PM2.5.FDMS),colour = 'red',linetype = 2)+
  geom_line(aes(y=ODIN.2.5f),colour = 'black',linetype = 1)+
  ylab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
tseries_pm2.5_1dy<-ggplot(data = merged.24hr, aes(x=date))+
  geom_line(aes(y=PM2.5.FDMS),colour = 'red',linetype = 2)+
  geom_line(aes(y=ODIN.2.5f),colour = 'black',linetype = 1)+
  ylab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol = 1)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
png('./corrected_odin_fdms_PM10.png',width = 2400,height = 1200)
tseries_pm10_1hr<-ggplot(data = all_merged.1hr, aes(x=date))+
  geom_line(aes(y=PM10.FDMS),colour = 'red',linetype = 2)+
  geom_line(aes(y=ODIN.10f),colour = 'black',linetype = 1)+
  ylab(bquote(PM[10]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
tseries_pm10_1dy<-ggplot(data = merged.24hr, aes(x=date))+
  geom_line(aes(y=PM10.FDMS),colour = 'red',linetype = 2)+
  geom_line(aes(y=ODIN.10f),colour = 'black',linetype = 1)+
  ylab(bquote(PM[10]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
grid.arrange(tseries_pm10_1hr,tseries_pm10_1dy,ncol = 1)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
png('./corrected_odin_fdms_PMcoarse.png',width = 2400,height = 1200)
tseries_pmcoarse_1hr<-ggplot(data = all_merged.1hr, aes(x=date))+
  geom_line(aes(y=PMcoarse),colour = 'red',linetype = 2)+
  geom_line(aes(y=ODIN.coarse),colour = 'black',linetype = 1)+
  ylab(bquote(PM[coarse]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
tseries_pmcoarse_1dy<-ggplot(data = merged.24hr, aes(x=date))+
  geom_line(aes(y=PMcoarse),colour = 'red',linetype = 2)+
  geom_line(aes(y=ODIN.coarse),colour = 'black',linetype = 1)+
  ylab(bquote(PM[coarse]~'('*mu~'g'~m^-3~')'))+
  theme(text=element_text(size=30))
grid.arrange(tseries_pmcoarse_1hr,tseries_pmcoarse_1dy,ncol = 1)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
png('./corrected_odin_fdms_PM.png',width = 1500,height = 750)
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,tseries_pmcoarse_1hr,tseries_pmcoarse_1dy,ncol=2)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
# System information
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: Fedora 21 (Twenty One)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] gridExtra_0.9.1 reshape2_1.4.1  openair_1.5     maps_2.3-9     
## [5] dplyr_0.4.1     lazyeval_0.1.10 ggplot2_1.0.1  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.11.6         mapdata_2.2-3       RColorBrewer_1.1-2 
##  [4] plyr_1.8.2          tools_3.2.0         digest_0.6.8       
##  [7] evaluate_0.7        gtable_0.1.2        nlme_3.1-120       
## [10] lattice_0.20-31     mgcv_1.8-6          png_0.1-7          
## [13] Matrix_1.2-0        DBI_0.3.1           mapproj_1.2-2      
## [16] yaml_2.1.13         parallel_3.2.0      hexbin_1.27.0      
## [19] proto_0.3-10        stringr_1.0.0       knitr_1.10.5       
## [22] cluster_2.0.1       RgoogleMaps_1.2.0.7 rmarkdown_0.6.1    
## [25] RJSONIO_1.3-0       latticeExtra_0.6-26 magrittr_1.5       
## [28] scales_0.2.4        htmltools_0.2.6     MASS_7.3-40        
## [31] assertthat_0.1      colorspace_1.2-6    labeling_0.3       
## [34] stringi_0.4-1       munsell_0.4.2
```


---
title: "odin_baseline_and_colocation_analysis.R"
author: "gustavo"
date: "Mon Jun 15 23:10:53 2015"
---
