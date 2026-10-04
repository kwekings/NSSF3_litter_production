
#############
# PHENOLOGY #
#############

library(mgcv)
library(MuMIn)
library(gamm4)

### Aggregate tree and plot basal area data ###
treeleaves <- lamina[which(!lamina$species %in% lianas),]
tl.sum <- aggregate(Dry.Mass ~ date1 + date2 + Plot, data=treeleaves, FUN=sum)
plot.type <- unique(lamina[c("Plot","Condition")])
tll <- merge(tl.sum, plot.type)
plot.ba <- with(tree, tapply(ba, plot, sum))
#pt <- c("D","W","D","D","W","W","D","D","W","W")
#boxplot(plot.ba ~ pt)
tll$ba <- plot.ba[match(tll$Plot, names(plot.ba))]

## add back twigs?
#tll <- merge(tll, twig, by.x=c("date1","Plot" ), by.y=c("date1","Plot"))
#tll$prod <- (tll$Dry.Mass.x + tll$Dry.Mass.y) / tll$ba*100^2 / 4
#tll$date2.y <- NULL; names(tll)[3] <- "date2"
## or don't add back twigs
tll$prod <- tll$Dry.Mass / tll$ba*100^2 / 4

summary(tll)
head(tll)

### Formatting weather data ###

location <- "D:\\National University of Singapore\\Chong Kwek Yan - CRSF\\Data\\Piezometer Data\\Khatib temperature\\"
pre <- "DAILYDATA_S122_"
months <- c("201806", "201807", gsub("-","", unique(tll$date2)))
weather <- list()
for(i in 1:length(unique(tll$date2))){
 weather[[i]] <- read.csv(paste0(location, pre, months[i],".csv"))}

library(stringr)
for(i in 1:length(weather)){
 weather[[i]]$Day <- str_pad(weather[[i]]$Day, 2, "left", "0")
 weather[[i]]$Month <- str_pad(weather[[i]]$Month, 2, "left", "0")
 weather[[i]]$date1 <- paste(weather[[i]]$Year, weather[[i]]$Month, weather[[i]]$Day, sep="-")}
weather <- do.call(rbind, weather)

names(weather)
weather <- weather[,c(9:11,14)]
names(weather) <- c("mean_temp", "max_temp", "min_temp","date1")

location <- "D:\\National University of Singapore\\Chong Kwek Yan - CRSF\\Data\\Piezometer Data\\Upp Peirce Res rainfall\\"
preUP <- "DAILYDATA_S69_"
preMD <- "DAILYDATA_S40_"
#months <- gsub("-","", unique(tll$date2))
rainfall <- list()
for(i in 1:length(unique(tll$date2))){
 if(months[i] %in% c("201903","201904","201905","201906","202008")){
  rainfall[[i]] <- read.csv(paste0(location, preMD, months[i],".csv"))
 } else {
  rainfall[[i]] <- read.csv(paste0(location, preUP, months[i],".csv"))
 } }

library(stringr)
for(i in 1:length(rainfall)){
 rainfall[[i]]$Day <- str_pad(rainfall[[i]]$Day, 2, "left", "0")
 rainfall[[i]]$Month <- str_pad(rainfall[[i]]$Month, 2, "left", "0")
 rainfall[[i]]$date1 <- paste(rainfall[[i]]$Year, rainfall[[i]]$Month, rainfall[[i]]$Day, sep="-")}
rainfall <- do.call(rbind, rainfall)
rainfall <- rainfall[,c(5,14)]
names(rainfall) <- c("total_rf","date1")

weather$total_rf <- rainfall$total_rf[match(weather$date1, rainfall$date1)]
weather$total_rf[which(weather$total_rf=="—")] <- 0
weather$total_rf <- as.numeric(weather$total_rf)
weather$date2 <- paste0(sapply(strsplit(weather$date1,"-"), "[[", 1),"-",sapply(strsplit(weather$date1,"-"), "[[", 2))
head(weather)
monthly_mean_temp <- with(weather, tapply(as.numeric(mean_temp), date2, mean, na.rm=T)) #NAs intro by coercion
weather$min_temp <- as.numeric(weather$min_temp) #NAs intro by coercion
weather$max_temp <- as.numeric(weather$max_temp) #NAs intro by coercion
index <- which(weather$mean_temp=="—")
for(i in index) weather$mean_temp[i] <- monthly_mean_temp[which(names(monthly_mean_temp)==weather$date2[i])]
weather$mean_temp <- as.numeric(weather$mean_temp)


### Importing weather data into tll data ###
weather$date1 <- as.Date(weather$date1)
tll$date1 <- as.Date(tll$date1)
# just use six vars: 
tll$rf <- rep(NA, nrow(tll))
tll$rf_lag <- rep(NA, nrow(tll))
tll$temp <- rep(NA, nrow(tll))
tll$temp_lag <- rep(NA, nrow(tll))
tll$temp_max <- rep(NA, nrow(tll))
tll$temp_min <- rep(NA, nrow(tll))

summary(weather)
for(i in 1:nrow(tll)){
 # rainfall of 30 days preceding collection date
 tll$rf[i] <- mean(weather$total_rf[which(
  weather$date1 < tll$date1[i] & weather$date1 > (tll$date1[i]-31)
 )])
 # rainfall of 30 days preceding one month before the collection date
 tll$rf_lag[i] <- mean(weather$total_rf[which(
  weather$date1 < (tll$date1[i]-30) & weather$date1 > (tll$date1[i]-61)
 )])
 # temp of 30 days preceding collection date
 tll$temp[i] <- mean(weather$mean_temp[which(
  weather$date1 < tll$date1[i] & weather$date1 > (tll$date1[i]-31)
 )])
 # temp of 30 days preceding one month before the collection date
 tll$temp_lag[i] <- mean(weather$mean_temp[which(
  weather$date1 < (tll$date1[i]-30) & weather$date1 > (tll$date1[i]-61)
 )])
 # min temp of the 30 days preceding collection date
 tll$temp_min[i] <- min(weather$min_temp[which(
  weather$date1 < tll$date1[i] & weather$date1 > (tll$date1[i]-31)
 )], na.rm=T)
 # max temp of the 30 days preceding collection date
 tll$temp_max[i] <- max(weather$max_temp[which(
  weather$date1 < tll$date1[i] & weather$date1 > (tll$date1[i]-31)
 )], na.rm=T)
}

summary(tll)
#pairs.cor(tll[,7:13])
#write.csv(tll, "data for GAMM.csv")

### START HERE ###
setwd("D:\\National University of Singapore\\Chong Kwek Yan - CRSF\\Past projects - DO NOT SHARE\\R_PJ\\publish")
tll <- read.csv("data for GAMM.csv", header=T)
tll$date1 <- as.Date(tll$date1)

### MODEL SELECTION ###
model.full<- gamm4(log(prod) ~ Condition + s(as.numeric(date1)) +
	rf + rf_lag + temp + temp_lag + 
	temp_max + temp_min +
	Condition:rf + Condition:rf_lag + Condition:temp + Condition:temp_lag + 
	Condition:temp_max + Condition:temp_min,
	random=~(1|Plot), data=tll)
summary(model.full$mer)
summary(model.full$gam)
plot(model.full$mer)
plot(model.full$gam)

uGamm.full <- uGamm(log(prod) ~ Condition + s(as.numeric(date1))
	+ rf + rf_lag + temp + temp_lag
	+ temp_max + temp_min
	+ Condition:rf + Condition:rf_lag + Condition:temp + Condition:temp_lag
	+ Condition:temp_max + Condition:temp_min
	, random=~(1|Plot), data=tll, lme4 = T)

dredged <- dredge(uGamm.full, m.lim=c(0,4), trace=T, rank="AICc",
	subset = !((rf && rf_lag) || (temp && temp_lag) ||
		(temp_max && temp_min) || (temp && temp_max) ||
		(temp && temp_min) || (temp_lag && temp_max) || (temp_lag && temp_min))
	)
dredged[dredged$delta<10,]
summary(get.models(dredged, subset=1)[[1]]$mer)
summary(get.models(dredged, subset=2)[[1]]$mer)
summary(get.models(dredged, subset=3)[[1]]$mer)

newdate <- sort(unique(tll$date1))

best <- gamm4(log(prod) ~ s(as.numeric(date1)) + Condition + temp,
	random=~(1|Plot), data=tll)
wpred <- predict(best$gam, newdata=data.frame(date1=newdate, Condition="W", temp=tll$temp[match(newdate,tll$date1)]), se.fit=T)
dpred <- predict(best$gam, newdata=data.frame(date1=newdate, Condition="D", temp=tll$temp[match(newdate,tll$date1)]), se.fit=T)

#plot(best$mer)
#plot(best$gam)
summary(best$gam)
summary(best$mer)

col.pal <- c("#9EC1A3" ,"#40798C",
	"#70A9A1", "#1F363D", "#8D8266")
col.t <- adjustcolor(col.pal, alpha.f=0.7)
#barplot(matrix(rep(1,5)), col=col.t)

tll2 <- tll[order(tll$date1),]

#pdf("Fig 2 ba-normalized revised.pdf", height=13, width=8, useDingbats=F)
#jpeg("Fig 2 ba-normalized new.jpg", height=13, width=8, units="in", res=300)
layout(matrix(1:2, ncol=1), heights=c(2,1))
par(mar=c(2,5.5,1.5,1), mgp=c(3,1,0))
plot(log(prod) ~ date1, data=tll, type="n",
	xlab="Date", cex.lab=1.5, cex=2, yaxt="n",
#	ylab=expression(paste("Basal area-normalized litter production (g", " fortnight"^" -1", " m"^" -4", ")")))
	ylab=expression(paste("Basal area-normalized litter production (g", " fortnight"^" -1", " m"^" -2", " m"^" -2", ")")))

polygon(c(wpred$fit+1.96*wpred$se.fit,
	rev(wpred$fit-1.96*wpred$se.fit)) ~ c(newdate,rev(newdate)), 
	border=F, col=col.t[2])
polygon(c(dpred$fit+1.96*dpred$se.fit,
	rev(dpred$fit-1.96*dpred$se.fit)) ~ c(newdate,rev(newdate)), 
	border=F, col=col.t[1])
lines(wpred$fit ~ newdate, lwd=8, col=col.pal[4], lty=2)
lines(dpred$fit ~ newdate, lwd=8, col=col.pal[5])
points(log(prod) ~ date1, data=tll, bg=ifelse(Condition=="W",col.t[4],col.t[5]), 
	col="white", lwd=2, cex=2.5, pch=ifelse(Condition=="W",21,23))
legend('topright', bty='n', cex=1.2, legend=c("Wet", "Dry"),
	title="Plot condition", pch=c(21,23), pt.cex=2, pt.lwd=2, 
	col="white", pt.bg=col.t[4:5])
legend('topright', bty='n', cex=1.2, legend=c("               ","            "),
	title="           ", lwd=4, lty=c(2,1), col=col.t[4:5])
axis(side=2, las=1, at=log(c(4,8,16,32,64,128)), labels=c(4,8,16,32,64,128))
mtext(side=3, adj=0, line=-1.2, cex=1.2, text=" a)")

par(mar=c(4.5, 5.5, 1, 1))
plot(temp ~ date1, data=tll2, lwd=3, type="l", las=1,
	xlim=range(tll$date1), ylim=c(26.5,29),
	ylab=expression(paste("Mean daily temperature (", degree, "C)")), xlab="Date", cex.lab=1.5)
mtext(side=3, adj=0, line=-1.2, cex=1.2, text=" b)")

dev.off()

r.squaredLR(best)
r.squaredLR(get.models(dredged, subset=1)[[1]])


summary(tll)
summary(lamina)





pairs.cor <- function (x,y,smooth=TRUE, digits=2,  ...)
{
  panel.cor <- function(x, y, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r.obj = cor.test(x, y,use="pairwise",...)
    r = as.numeric(r.obj$estimate)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(txt)
    cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex=cex*abs(r))
  }
panel.hist <- function(x)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="grey")
  }
pairs(x,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.smooth, ...)
}













