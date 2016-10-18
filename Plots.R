spVL.run <- read.table("ViralLoadOutput.txt", header=TRUE)
patient.run <- read.table("PatientRegistryOutput.txt", header=TRUE)

#tfinal <-365*30

#xgraphlimit = round(round(10*((tfinal+365)/365), digits = 0)/10)

# edits to patient.run file
patient.run <- subset(patient.run, patient.run$TimeInf > 0)

par(mfrow=c(3,3))


#######  1.  plot of relationship between VIRAL LOAD AND TRANSMISSION 

patient.run <- transform(patient.run, DonorsVLatTrans=log10(DonorsVLatTrans))
hist(patient.run$DonorsVLatTrans, breaks=50, freq=FALSE, las=1, xlab="VL at transmission (log10)", ylab="Frequency", main="VL at transmission", xlim=c(1,8))



#######  2.  plot of TIME OF INFECTION AND NUMBER OF TRANSMISSIONS
 
t_acute <- 91
 
new.list.1 <- patient.run[, "Donors_Total_Time_Inf_At_Trans"]
new.list.1 <- as.data.frame(new.list.1)
new.list.2 <- subset(new.list.1, new.list.1 > 0)
new.list.2 <- as.numeric(new.list.2$new.list.1)
new.list.2 <- as.data.frame(new.list.2)
new.list.acute <- subset(new.list.2, new.list.2 <= t_acute)
new.list.chronic <- subset(new.list.2, new.list.2 > t_acute)
length(new.list.acute$new.list.2) -> acute
length(new.list.chronic$new.list.2) -> not.acute
acute/not.acute -> freq.acute
hist(new.list.2$new.list.2/365, breaks=60.3, freq=FALSE, xlab="Years after infection", ylab="Density", main="Time of infection at transmission", 
     xlim=c(0,15), las=1)
text(5, 0.15, paste("Freq inf in acute infection =", format(freq.acute[1], digits=2)), adj=c(0,0), col=2, cex=0.9)



########  3.  plot of HERITABILITY OF SET POINT VIRAL LOAD
 
test <- subset(patient.run, select=c("TimeInf", "NumRecipients", "LogSP0", "Donors_LogSP0", "Total_Time_Inf"))

x.365 <- numeric(0)

for (i in 0:(xgraphlimit-2)) {
patient.run.yr <- patient.run[patient.run$TimeInf > i*365 & patient.run$TimeInf < (i+1)*365, ]
h <- (lm(patient.run.yr$LogSP0 ~ patient.run.yr$Donors_LogSP0))
d <- format(h$coefficients[2], digits=2)
x.365 <- c(x.365, d)
}

plot

years <- c(1:(xgraphlimit-1))
test = data.frame(v1=x.365, years=years)
test$v1 <- as.numeric(as.character(test$v1))

plot(test$years, test$v1, type="o", pch=18, lty=1, col="red", bg="red", xlim=c(0, tfinal/365), ylim=c(0, 1.0), xlab="Years", 
ylab="Regression coefficient", las=1, main="Heritability")

  
  
########  4.  plot of Time To AIDS ~ SPVL

aa <- subset(patient.run, patient.run$Status=="DiedAIDS" | patient.run$Status=="Alive(HIV+)")

plot(aa$LogSP0, aa$TimeToAIDS/365, xlab="SPVL", ylab="Time to AIDS, years", main="SPVL and time to AIDS", las=1, xlim=c(2,8), ylim=c(0, 40))

lines(lowess(aa$LogSP0, aa$TimeToAIDS/365, f=0.01), col="blue", lwd=2)

abline(10, 0, lty=2)
test.a <- subset(aa, aa$LogSP0 > 2 & aa$LogSP0 < 3)
x <- mean(test.a$TimeToAIDS, na.rm=TRUE)/365
points(2.5, x, col="red")

abline(10, 0, lty=2)


   
  
#######  5.  plot of DISTRIBUTIONS OF SET POINT VIRAL LOAD FOR ALL INDIVIDUALS AT CERTAIN YEARS OF THE EPIDEMIC

patient.run <- read.table("PatientRegistryOutput.txt", header=TRUE)
 
spvl.0 <- subset(patient.run, patient.run$TimeInf/365 < 0)
tstopped = max(spVL.run$time)
if (tstopped/365 > 9) spvl.10 <- subset(patient.run, patient.run$TimeInf/365 > 9 & patient.run$TimeInf/365 < 10)
if (tstopped/365 > 19) spvl.20 <- subset(patient.run, patient.run$TimeInf/365 > 19 & patient.run$TimeInf/365 < 20)

x.0 <- spvl.0$LogSP0
if (tstopped/365 > 9) x.10 <- spvl.10$LogSP0
if (tstopped/365 > 19) x.20 <- spvl.20$LogSP0

x.0.density <- density(x.0, bw=0.4)
if (tstopped/365 > 9) x.10.density <- density(x.10, bw=0.3)
if (tstopped/365 > 19) x.20.density <- density(x.20, bw=0.3)

plot(x.0.density, col="black", lwd=2, xlab="SPVL", ylab="Density", main="SPVL distributions", xlim = c(2,8), 
ylim=c(0,0.6), las=1)
if (tstopped/365 > 9) lines(x.10.density, col="blue", lwd=2)
if (tstopped/365 > 19) lines(x.20.density, col="orange", lwd=2)

text(2, .5, c("time zero"), cex=0.9, col="black", pos=4)
text(2, .4, c("year 10"), cex=0.9, col="blue", pos=4)
text(2, .3, c("year 20"), cex=0.9, col="orange", pos=4)

   
  
#######  6.  plot of SET POINT VIRAL LOAD OVER TIME
 
patient.run <- subset(patient.run, patient.run$TimeInf > 0)
 
spVL.time.inf.years <- patient.run$TimeInf/365

plot(spVL.time.inf.years, patient.run$LogSP0[patient.run$LogSP0 != 0], col="white", xlab="Years", ylab="Average SPVL",
    main="Average SPVL through time", las=1, xlim=c(0, xgraphlimit), ylim=c(3, 5))

c <- smooth.spline(spVL.time.inf.years, patient.run$LogSP0)

abline(h=seq(3, 5, by=0.5), col="light gray", ylim=c(0, tfinal/365))

#par(new=TRUE)

lines(lowess(spVL.time.inf.years, patient.run$LogSP0, f=1), lty = 1, lwd=2, col = "red")
#lines(lowess(spVL.run$time[spVL.run$time > 0]/365, spVL.run$AveLogV[spVL.run$time > 0], f=0.1), lty=1, lwd=2, col="blue")
 
#text(10, 4.4, c("Mean SPVL"), cex=0.9, col="red", pos=4)
#text(15, 3.5, c("(Mean) CVL"), cex=0.9, col="blue", pos=4)

 


#######  7.  plot of PREVALENCE OVER TIME
 
prev <- spVL.run$Infected/(spVL.run$Susc + spVL.run$Infected) 
plot(spVL.run$time/365, prev, col="red", 
	xlab="Years", ylab="Prevalence", main="Prevalence", xlim=c(0, xgraphlimit), 
	ylim=c(0, max(prev)), las=1, lwd=2)



#######  8.  plot of INCIDENCE OVER TIME
 
#install.packages("zoo")
library(zoo)
#install.packages("plyr")
library(plyr)

#aggregate times of infection by year
#make a column of X years
#aggregate by year

### Incidence (new infections per year)

test.spvl <- subset(spVL.run, select=c("time", "Susc", "Infected", "dAIDS", "AveLogV", "Vstd"))
test.spvl <- test.spvl[test.spvl$time > 0, ]
test.spvl["year"] <- NA
test.spvl$year <- ceiling((test.spvl$time)/365)

test.patient <- subset(patient.run, select=c("TimeInf", "LogSP0"))
test.patient <- test.patient[test.patient$TimeInf > 0, ]
test.patient["year"] <- NA
test.patient <- test.patient[order(test.patient$TimeInf) , ]
test.patient$year <- ceiling((test.patient$TimeInf)/365)

#Aggregate mean CVLs by year
#Count times of infection by year (just counting the number of rows for each year)
#aggregation steps:

mean.cvl.by.year <- aggregate(test.spvl, by=list(test.spvl$year), FUN=mean)
mean.cvl.by.year <- mean.cvl.by.year[-nrow(mean.cvl.by.year), ]

inf.by.year <- count(test.patient, "year")
inf.by.year <- inf.by.year[-nrow(inf.by.year), ]

susc.by.year <- mean.cvl.by.year$Susc
inc.by.year <- inf.by.year$freq/susc.by.year
tenyears <- cbind(susc.by.year, inf.by.year, inc.by.year)
tenyears <- subset(tenyears, select=c("year", "inc.by.year"))

plot(tenyears$year, 100*tenyears$inc.by.year, las=1, col="red", main="Incidence through time", 
	ylab="Incidence, per 100 person years", xlab="Years", xlim=c(0, xgraphlimit), 
	ylim=c(0, max(100*tenyears$inc.by.year)))

