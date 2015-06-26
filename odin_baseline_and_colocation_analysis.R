#' ---
#' title: "ODIN baseline and colocation analysis - June 2015"
#' author: Gustavo Olivares
#' output:
#'  html_document:
#'    keep_md: true
#' ---

###############################################################################
# Analysis script for the manuscript "The Outdoor Dust Information Node (ODIN)
# - Development and performance assessment of a low cost ambient dust sensor" 
# submitted to Atmospheric Measurement Techniques. June 2015
###############################################################################
## Preamble ####
# Libraries
require(ggplot2)
require(openair)
require(reshape2)
require(gridExtra)
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
dev.off()
long_scatter_temperature_plot
summary(Temp.lm<-lm(data=long_scatter_temperature,Dust~Temperature))
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
all_merged.1hr$ODIN_drift[length(all_merged.1hr$ODIN_drift)] # Baseline at the end 

# Remove the baseline drift from the raw ODIN data
all_merged.1hr$ODIN_raw <- all_merged.1hr$ODIN
all_merged.1hr$ODIN_detrend<-all_merged.1hr$ODIN_raw - all_merged.1hr$ODIN_drift
# Temperature correction
summary(odin_T<-lm(data=all_merged.1hr,ODIN_detrend~Temperature,subset = Temperature>10))

format(coef(odin_T),digits = 2)
format((confint(odin_T)[,2]-confint(odin_T)[,1])/2,digits = 1)

png('./odin_T_scatter.png',width = 750,height = 750)
odin_T_scatter<-qplot(Temperature, ODIN_detrend,data = all_merged.1hr)+
  #geom_point(colour = 'blue')
  geom_abline(intercept = coef(odin_T)[1],slope = coef(odin_T)[2])+
  ylab(bquote(ODIN[detrended]~'(mV)'))+
  xlab('Temperature (C)')+
  theme_bw()+
  theme(text=element_text(size=30))
odin_T_scatter
dev.off()
odin_T_scatter
all_merged.1hr$ODIN <- all_merged.1hr$ODIN_detrend - predict(odin_T,newdata = all_merged.1hr)

# Dust performance using ECan data for calibration
# Calibration expression:
#  $ODIN_{calibrated}=A*ODIN_{raw}+B$

# Using only the first 30% of the data to fit the models
data_selection.1hr <- rbinom(length(all_merged.1hr$date),1,0.3)
data_selection.1hr <- (seq(length(all_merged.1hr$date))<(length(all_merged.1hr$date)/3))
# PM$_{2.5}$ fdms
summary(odin1.lm.1hr.pm2.5.fdms<-
          lm(data=all_merged.1hr,PM2.5.FDMS~ODIN,
             subset = data_selection.1hr == 1))
format(coef(odin1.lm.1hr.pm2.5.fdms),digits = 1)
format((confint(odin1.lm.1hr.pm2.5.fdms)[,2]-confint(odin1.lm.1hr.pm2.5.fdms)[,1])/2,digits = 1)
### Calculate the corrected data
all_merged.1hr$ODIN.2.5f<-predict(odin1.lm.1hr.pm2.5.fdms,
                                  newdata = all_merged.1hr)
sqrt( mean( (all_merged.1hr$ODIN.2.5f-all_merged.1hr$PM2.5.FDMS)^2 , na.rm = TRUE ) ) #RMSE for the whole dataset
cor(all_merged.1hr$ODIN.2.5f,all_merged.1hr$PM2.5.FDMS,use = "pairwise")^2 # R^2 for the whole dataset
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
dev.off()
grid.arrange(tseries_pm2.5,tseries_ODIN,ncol = 1)

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
dev.off()
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol = 1)

png('./corrected_odin_fdms_PM.png',width = 1500,height = 1500)
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol=2)
dev.off()
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,ncol=2)

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
dev.off()
scatter_fitted_PM2.5_hour
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
scatter_fitted_PM2.5_day

png('./tseries_and_scatter.png',width = 1500, height = 1500)
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,scatter_fitted_PM2.5_hour,scatter_fitted_PM2.5_day,ncol=2)
dev.off()
grid.arrange(tseries_pm2.5_1hr,tseries_pm2.5_1dy,scatter_fitted_PM2.5_hour,scatter_fitted_PM2.5_day,ncol=2)


# System information
sessionInfo()









