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
require(dplyr)
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
long_scatter_temperature_plot<-ggplot(long_scatter_temperature, aes(x=Temperature, y=Dust, color=Sensor)) +
  geom_point(shape=19,alpha=1/3,size = 10)+
  theme_bw()+
  theme(text = element_text(size = 30))+
  xlab('Temperature (C)')+
  ylab('Sensor response (mV)')
long_scatter_temperature_plot
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
long_scatter_temperature_plot
```

```
## Warning in loop_apply(n, do.ply): Removed 744 rows containing missing
## values (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-1.png) 

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

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-2.png) 

```r
odin01_lag=lag_test$lag[which.max(lag_test$acf)]
odin_01$date=odin_01$date-odin01_lag*10*60
odin_01$ODIN = odin_01$Dust
odin_01$Dust<-NULL
all_merged_FULL<-merge(odin_01,ecan_data,by = 'date',all = TRUE)
all_merged.10min<-timeAverage(all_merged_FULL,avg.time = '10 min')
# Full dataset 1 hour  
all_merged.1hr<-timeAverage(all_merged.10min,avg.time='1 hour')
all_merged.1hr$PMcoarse<-all_merged.1hr$PM10.FDMS - all_merged.1hr$PM2.5.FDMS

# Remove drift from ODIN raw data
# Estimate the baseline from a simple linear regression
all_merged.1hr$ODIN_drift<-predict(lm(all_merged.1hr$ODIN~seq(all_merged.1hr$ODIN)),newdata = all_merged.1hr)

all_merged.1hr$ODIN_drift[1] # Baseline at the beginning
```

```
## [1] 403.3056
```

```r
all_merged.1hr$ODIN_drift[length(all_merged.1hr$ODIN_drift)] # Baseline at the end 
```

```
## [1] 287.552
```

```r
# Remove the baseline drift from the raw ODIN data
all_merged.1hr$ODIN_raw <- all_merged.1hr$ODIN
all_merged.1hr$ODIN_detrend<-all_merged.1hr$ODIN_raw - all_merged.1hr$ODIN_drift
# Temperature correction
summary(odin_T<-lm(data=all_merged.1hr,ODIN_detrend~Temperature,subset = Temperature>10))
```

```
## 
## Call:
## lm(formula = ODIN_detrend ~ Temperature, data = all_merged.1hr, 
##     subset = Temperature > 10)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -62.750  -7.120  -0.781   9.946  28.070 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -49.2798     3.8994  -12.64   <2e-16 ***
## Temperature   4.3821     0.2553   17.16   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13.49 on 189 degrees of freedom
##   (22 observations deleted due to missingness)
## Multiple R-squared:  0.6092,	Adjusted R-squared:  0.6071 
## F-statistic: 294.6 on 1 and 189 DF,  p-value: < 2.2e-16
```

```r
format(coef(odin_T),digits = 2)
```

```
## (Intercept) Temperature 
##     "-49.3"     "  4.4"
```

```r
format((confint(odin_T)[,2]-confint(odin_T)[,1])/2,digits = 1)
```

```
## (Intercept) Temperature 
##       "7.7"       "0.5"
```

```r
png('./odin_T_scatter.png',width = 750,height = 750)
odin_T_scatter<-qplot(Temperature, ODIN_detrend,data = all_merged.1hr)+
  #geom_point(colour = 'blue')
  geom_abline(intercept = coef(odin_T)[1],slope = coef(odin_T)[2])+
  ylab(bquote(ODIN[detrended]~'(mV)'))+
  xlab('Temperature (C)')+
  theme_bw()+
  theme(text=element_text(size=30))
odin_T_scatter
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_point).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
odin_T_scatter
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-3.png) 

```r
all_merged.1hr$ODIN <- all_merged.1hr$ODIN_detrend - predict(odin_T,newdata = all_merged.1hr)

# Dust performance using ECan data for calibration
# Calibration expression:
#  $ODIN_{calibrated}=A*ODIN_{raw}+B$

# Using only the first 1/3 of the data to fit the models
data_selection.1hr <- (seq(length(all_merged.1hr$date))<(length(all_merged.1hr$date)/3))
# PM$_{2.5}$ fdms
summary(odin1.lm.1hr.pm2.5.fdms<-
          lm(data=all_merged.1hr,PM2.5.FDMS~ODIN,
             subset = data_selection.1hr == 1))
```

```
## 
## Call:
## lm(formula = PM2.5.FDMS ~ ODIN, data = all_merged.1hr, subset = data_selection.1hr == 
##     1)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -34.692  -6.487  -0.585   4.417  41.715 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 18.73469    0.89114   21.02   <2e-16 ***
## ODIN         0.63465    0.02067   30.70   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 10.95 on 158 degrees of freedom
##   (15 observations deleted due to missingness)
## Multiple R-squared:  0.8564,	Adjusted R-squared:  0.8555 
## F-statistic: 942.6 on 1 and 158 DF,  p-value: < 2.2e-16
```

```r
format(coef(odin1.lm.1hr.pm2.5.fdms),digits = 1)
```

```
## (Intercept)        ODIN 
##      "18.7"      " 0.6"
```

```r
format((confint(odin1.lm.1hr.pm2.5.fdms)[,2]-confint(odin1.lm.1hr.pm2.5.fdms)[,1])/2,digits = 1)
```

```
## (Intercept)        ODIN 
##      "1.76"      "0.04"
```

```r
### Calculate the corrected data
all_merged.1hr$ODIN.2.5f<-predict(odin1.lm.1hr.pm2.5.fdms,
                                  newdata = all_merged.1hr)
sqrt( mean( (all_merged.1hr$ODIN.2.5f-all_merged.1hr$PM2.5.FDMS)^2 , na.rm = TRUE ) ) #RMSE for the whole dataset
```

```
## [1] 15.41
```

```r
cor(all_merged.1hr$ODIN.2.5f,all_merged.1hr$PM2.5.FDMS,use = "pairwise")^2 # R^2 for the whole dataset
```

```
## [1] 0.7242653
```

```r
# 24 hour dataset
merged.24hr <- timeAverage(all_merged.1hr,avg.time = '24 hour',statistic = 'mean')

png('./raw_odin_fdms.png',width = 1500,height = 1500)
tseries_pm2.5<-ggplot(data = all_merged.1hr, aes(x=date,y=PM2.5.FDMS))+
  geom_line(colour = 'red')+
  ylab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))+
  theme_bw()+
  theme(text=element_text(size=30))
tseries_ODIN<-ggplot(data = all_merged.1hr, aes(x=date,y=ODIN_raw))+
  geom_line(colour = 'black')+
  theme_bw()+
  ylab('ODIN(mV)')+
  theme(text=element_text(size=30))
grid.arrange(tseries_pm2.5,tseries_ODIN,ncol = 1)
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
grid.arrange(tseries_pm2.5,tseries_ODIN,ncol = 1)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-4.png) 

```r
png('./corrected_odin_fdms_PM2.5.png',width = 2400,height = 1200)
tseries_pm2.5_1hr<-ggplot(data = all_merged.1hr, aes(x=date))+
  geom_line(aes(y=PM2.5.FDMS,colour = '0'))+
  geom_line(aes(y=ODIN.2.5f,colour = '1'))+
  ylab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))+
  theme_bw()+
  ggtitle('a)')+
  scale_colour_manual(values = c('red','black'),name='',labels = c(expression(PM[2.5]),'ODIN'))+
  theme(text=element_text(size=30),legend.position=c(0.5,0.9))
tseries_pm2.5_1dy<-ggplot(data = merged.24hr, aes(x=date))+
  geom_line(aes(y=PM2.5.FDMS,colour = '0'))+
  geom_line(aes(y=ODIN.2.5f,colour = '1'))+
  ylab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))+
  theme_bw()+
  ggtitle('b)')+
  scale_colour_manual(values = c('red','black'),name='',labels = c(expression(PM[2.5]),'ODIN'))+
  theme(text=element_text(size=30),legend.position=c(0.5,0.9))
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
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol = 1)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-5.png) 

```r
png('./corrected_odin_fdms_PM.png',width = 1500,height = 1500)
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol=2)
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
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol=2)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-6.png) 

```r
png('./corrected_odin_scatter_hr.png',width = 1500,height = 1500)
scatter_fitted_PM2.5_hour<-ggplot(data = all_merged.1hr,aes(x=PM2.5.FDMS,y = ODIN.2.5f))+
  geom_point(shape = 1)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  ggtitle('c)')+
  theme(text=element_text(size=30))+
  ylab(bquote(ODIN~estimate~'('*mu~'g'~m^-3~')'))+
  xlab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))
scatter_fitted_PM2.5_hour
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_point).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
scatter_fitted_PM2.5_hour
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-7.png) 

```r
png('./corrected_odin_scatter_dy.png',width = 1500,height = 1500)
scatter_fitted_PM2.5_day<-ggplot(data = merged.24hr,aes(x=PM2.5.FDMS,y = ODIN.2.5f))+
  geom_point(shape = 1)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  ggtitle('d)')+
  theme(text=element_text(size=30))+
  ylab(bquote(ODIN~estimate~'('*mu~'g'~m^-3~')'))+
  xlab(bquote(PM[2.5]~'('*mu~'g'~m^-3~')'))
scatter_fitted_PM2.5_day
dev.off()
```

```
## png 
##   2
```

```r
scatter_fitted_PM2.5_day
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-8.png) 

```r
png('./tseries_and_scatter.png',width = 1500, height = 1500)
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,scatter_fitted_PM2.5_hour,scatter_fitted_PM2.5_day,ncol=2)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_point).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,scatter_fitted_PM2.5_hour,scatter_fitted_PM2.5_day,ncol=2)
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_path).
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing missing values
## (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-9.png) 

```r
### Inlet performance
# Load wind speed and direction data

wind_data<-read.table("http://files.figshare.com/2166068/StAlbans24July_14August_2014_WIND.tsv",
                      header = TRUE,
                      sep="\t",
                      quote = "\"'",
                      dec=".")
names(wind_data)<-c('date','wd','ws')
wind_data$date<-as.POSIXct(wind_data$date,format = "%d/%m/%Y %I:%M:%S %p")
wind_merge <- merge(all_merged.1hr,wind_data, by = 'date', all=TRUE)
wind_merge <- data.frame(wind_merge,bin_ws=cut(wind_merge$ws,c(0,1,2,3,4,10),include.lowest = TRUE))
wind_merge$ERROR <- wind_merge$PM2.5.FDMS - wind_merge$ODIN.2.5f

# Wind direction
inlet_wd<-ggplot(wind_merge,aes(x = wd,y = ERROR)) + 
  geom_point(shape=19,alpha=1/3,size = 10) +
  theme_bw()+
  ggtitle('a)')+
  theme(text=element_text(size=30))+
  ylab(bquote(PM[2.5]-ODIN~'('*mu~'g'~m^-3~')'))+
  xlab('Wind Direction') +
  scale_x_continuous(breaks = c(0,45,90,135,180,225,270,315),
                     labels = c('N','NE','E','SE','S','SW','W','NW'))
inlet_wd
```

```
## Warning in loop_apply(n, do.ply): Removed 23 rows containing missing values
## (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-10.png) 

```r
# Wind speed
inlet_ws<-ggplot(wind_merge,aes(x = ws,y = ERROR)) + 
  geom_point(shape=19,alpha=1/3,size = 10) +
  theme_bw()+
  ggtitle('b)') +
  theme(text=element_text(size=30)) +
  ylab(bquote(PM[2.5]-ODIN~'('*mu~'g'~m^-3~')')) +
  xlab(bquote(Wind~Speed~'(m '~s^-1~')')) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))
inlet_ws
```

```
## Warning in loop_apply(n, do.ply): Removed 23 rows containing missing values
## (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-11.png) 

```r
png('./inlet_scatter.png',width = 1500, height = 1500)
grid.arrange(inlet_wd,inlet_ws,ncol=1)
```

```
## Warning in loop_apply(n, do.ply): Removed 23 rows containing missing values
## (geom_point).
```

```
## Warning in loop_apply(n, do.ply): Removed 23 rows containing missing values
## (geom_point).
```

```r
dev.off()
```

```
## png 
##   2
```

```r
grid.arrange(inlet_wd,inlet_ws,ncol=1)
```

```
## Warning in loop_apply(n, do.ply): Removed 23 rows containing missing values
## (geom_point).
```

```
## Warning in loop_apply(n, do.ply): Removed 23 rows containing missing values
## (geom_point).
```

![](odin_baseline_and_colocation_analysis_files/figure-html/unnamed-chunk-1-12.png) 

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
date: "Wed Jul  8 12:05:44 2015"
---
