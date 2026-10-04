#########################
# LOAD DATA & SETUP DFs #
#########################

setwd("D:\\National University of Singapore\\Chong Kwek Yan - CRSF\\Past projects - DO NOT SHARE\\R_PJ\\publish")

# Lamina
lamina <- read.csv("Cleaned lamina+PJ.csv", header = T, stringsAsFactors=T)
lamina$date1 <- as.POSIXct(lamina$General.date, format ="%d/%m/%Y", tz="GMT")
lamina$date2 <- format(as.Date(lamina$date1), "%Y-%m")
summary(lamina)
lianas <- levels(lamina$species)[c(7,9,15,23,24,30,33,36,37,39,57,60,108,109,111)]

# Twigs
twig <- read.csv("Cleaned twigs+PJ.csv", header=T)
summary(twig)
twig$date1 <- as.POSIXct(twig$General.date, format ="%d/%m/%Y", tz="GMT")
twig$date2 <- format(as.Date(twig$date1), "%Y-%m")
summary(twig)

####################
# LITTER COMM NMDS #
####################

spp.by.plot <- with(lamina, tapply(Dry.Mass, list(Plot, species), sum))
spp.by.plot[is.na(spp.by.plot)] <- 0

library(vegan)

leaf.nmds <- metaMDS(spp.by.plot, dist="bray", k=2)
plot(leaf.nmds, type="t")

wetplots<-c("Q10","Q3","Q4","Q6","Q9")

# select only the IDed trees from this
rownames(leaf.nmds$species)
sel <- c(6,8,10:14,16:19,21,22,25:28,31,32,34,35,38,40:46,49:55,58,59,61,110,112)
hist(cex.leaf <- log(apply(spp.by.plot,2,sum)[sel]))

spp.sel <- rownames(leaf.nmds$species)[sel]
identical(spp.sel, colnames(spp.by.plot)[sel])

#jpeg("Leaf litter communities.jpg", width=14, height=8, units="in", res=300)
par(mar=c(5,5,2,2))
plot(leaf.nmds$species, ylab="NMDS 2", xlab="NMDS 1", cex.lab=2, cex.axis=1.5, type="n", xlim=c(-1.5,1), yaxt="n")
axis(side=2, at=seq(-0.8,0.4,0.4), cex.axis=1.5, las=1)
text(leaf.nmds$species[sel,], labels=rownames(leaf.nmds$species[sel,]), cex=cex.leaf/4, col="grey50", font=3)
points(leaf.nmds$points, pch=16, cex=7, col=ifelse(rownames(leaf.nmds$points) %in% wetplots, "#41658AB3", "#70A37FB3"))
text(leaf.nmds$points, labels=rownames(leaf.nmds$points), col="white")
legend('bottomleft', bty='n', cex=1.5, legend=c("Wet", "Dry"),
	title="Plot condition", pch=16, pt.cex=3,	col=c("#41658AB3","#70A37FB3"))
dev.off()

# PERMANOVA
pt <- factor(ifelse(rownames(leaf.nmds$points) %in% wetplots, "wet", "dry"))
adonis(spp.by.plot ~ pt, dist="bray", permutations=99999)

##################
# TREE COMM NMDS #
##################

tree <- read.csv("ten plot trees.csv", header=T)
summary(tree)
tree <- droplevels(na.omit(tree))
tree <- droplevels(tree[tree$species != "",])
tree$ba <- pi*(tree$dbh_2018/2)^2
summary(tree)


library(reshape2)
tree.comm <- dcast(data=tree, formula=plot~species, fun.aggregate=sum, value.var="ba", fill=0)
rownames(tree.comm) <- tree.comm$plot
tree.comm$plot <- NULL

tree.nmds <- metaMDS(tree.comm, k=2, dist="bray")
#plot(tree.nmds, type="t")

spp.centroids <- tree.nmds$species[which(apply(tree.comm, 2, sum) > 1000),]
spp.cex <- apply(tree.comm, 2, sum)[which(apply(tree.comm, 2, sum) > 1000)]
spp.cex <- sqrt(spp.cex)/40

abbrev <- paste0(
  substr(sapply(strsplit(rownames(spp.centroids), " "), "[[", 1),0,3),
  ".",
  substr(sapply(strsplit(rownames(spp.centroids), " "), "[[", 2),0,3) )

#jpeg("tree communities abbrev.jpg", width=14, height=8, units="in", res=300)
par(mar=c(5,5,2,2))
plot(tree.nmds$species, ylab="NMDS 2", xlab="NMDS 1", cex.lab=2, cex.axis=1.5, type="n", xlim=c(-0.6,0.8), ylim=c(0.65,-0.4))
points(tree.nmds$points, cex=7, 
  col=ifelse(rownames(tree.nmds$points) %in% wetplots, "#41658AB3", "#70A37FB3"),
  pch=ifelse(rownames(tree.nmds$points) %in% wetplots, 16, 17))
text(jitter(spp.centroids[which(spp.cex>0.95),],250), labels=abbrev[which(spp.cex>0.95)], cex=spp.cex[which(spp.cex>0.95)], col="grey50", font=3)
#text(tree.nmds$points, labels=rownames(tree.nmds$points), col="white")
legend('bottomleft', bty='n', cex=1.5, legend=c("Swamp", "Non-swamp"),
	title="Plot type", pch=c(16,17), pt.cex=3,	col=c("#41658AB3","#70A37FB3"))
dev.off()

# combined figure

#jpeg("combined ordinations abbrev.jpg", width=10, height=12, units="in", res=300)
#pdf("combined ordinations abbrev.pdf", width=10, height=12, useDingbats=F)
par(mfrow=c(2,1), oma=c(3,3,0,0), mar=c(2,3,2,2))

plot(tree.nmds$species, ylab="NMDS 2", xlab="NMDS 1", cex.lab=2, cex.axis=1.5, type="n", xlim=c(-0.6,0.7), ylim=c(0.6,-0.4))
text(spp.centroids, labels=abbrev, cex=spp.cex, col="grey50", font=3)
points(tree.nmds$points, pch=16, cex=7, col=ifelse(rownames(tree.nmds$points) %in% wetplots, "#41658AB3", "#70A37FB3"))
text(tree.nmds$points, labels=rownames(tree.nmds$points), col="white")
mtext(side=3, line=-2, text=" a) Canopy trees", cex=2, adj=0)

plot(leaf.nmds$species, ylab="", xlab="", cex.lab=2, cex.axis=1.5, type="n", xlim=c(-1.2,1), yaxt="n")
axis(side=2, at=seq(-0.8,0.4,0.4), cex.axis=1.5, las=1)
text(leaf.nmds$species[sel,], labels=spp.sel, cex=cex.leaf/4, col="grey50", font=3)
points(leaf.nmds$points, pch=16, cex=7, col=ifelse(rownames(leaf.nmds$points) %in% wetplots, "#41658AB3", "#70A37FB3"))
text(leaf.nmds$points, labels=rownames(leaf.nmds$points), col="white")
mtext(side=3, line=-2, text=" b) Leaf litter", cex=2, adj=0)

legend('bottomleft', bty='n', cex=1.5, legend=c("Wet", "Dry"),
	title="Plot condition", pch=16, pt.cex=4,	col=c("#41658AB3","#70A37FB3"))

mtext(side=1, outer=T, text="NMDS 1", cex=2, line=1)
mtext(side=2, outer=T, text="NMDS 2", cex=2, line=1)

dev.off()


### summary stats
length(tree.comm)
nrow(leaf.nmds$species)

#######################################
# SPECIES INTRINSIC LITTER PRODUCTION #
#######################################

spp.sel
# 5 species (not present in 10 plots) omitted in this step (n=36):
ba.by.plot <- tree.comm[,na.omit(match(spp.sel,colnames(tree.comm)))]

# extract litter collection from these 36 species
ind.lam <- which(colnames(spp.by.plot) %in% colnames(ba.by.plot))
litter.consol <- melt(spp.by.plot[,ind.lam], varnames=c("plot", "species"), value.name = "litter.collection")

# combine with basal area 
tpba <- melt(ba.by.plot, value.name = "ten.ba")
litter.consol$ten.ba <- tpba$ten.ba

# survey duration
for(i in 1:length(spp.sel)){
	litter.consol$start[which(litter.consol$species==spp.sel[i])] <- 
	  as.character(min(lamina[which(lamina$species==spp.sel[i]),"date1"], na.rm=T))
	if(max(lamina[which(lamina$species==spp.sel[i]),"date1"], na.rm=T) ==
	  min(lamina[which(lamina$species==spp.sel[i]),"date1"], na.rm=T) ){
		litter.consol$duration[which(litter.consol$species==spp.sel[i])] <- 
		  as.numeric(max(lamina[,"date1"], na.rm=T) - 
		  min(lamina[which(lamina$species==spp.sel[i]),"date1"], na.rm=T))
	  } else {		
	litter.consol$duration[which(litter.consol$species==spp.sel[i])] <- 
		as.numeric(max(lamina[which(lamina$species==spp.sel[i]),"date1"], na.rm=T) - 
		  min(lamina[which(lamina$species==spp.sel[i]),"date1"], na.rm=T))
	  }
}

# amount of litter produced per year
litter.consol$litter.production <- with(litter.consol, 
	litter.collection/4 	# convert units to g/m2 of plot area
	/ (duration/365))		# per year of survey

# abbreviate species name
spp.char <- strsplit(as.character(litter.consol$species), " ")
litter.consol$code <- 
  toupper(paste0(
	substr(sapply(spp.char, "[[", 1), 0, 1),
	substr(sapply(spp.char, "[[", 2), 0, 2)
  ))

#hist(litter.consol$litter.production)

plot(log(litter.production+1) ~ log(ten.ba), data=litter.consol, type="n")
text(log(litter.production+1) ~ log(ten.ba), data=litter.consol, labels=code,
  col=rainbow(36)[as.numeric(litter.consol$species)], cex=1.5)

summary(litter.consol)
unique(litter.consol$species)
# remove plots in which no litter was collected, or where species did not occur
litter.consol <- litter.consol[-which(litter.consol$litter.production==0),]
# the next line removes 27 data points and 2 species: Elaeocarpus stipularis and Callophyllum wallichianum
litter.consol <- litter.consol[-which(litter.consol$ten.ba==0),]
summary(litter.consol)
#hist(log(litter.consol$litter.production))
#hist(log(litter.consol$ten.ba))

##########################
# wet-dry specialization #
##########################

SSI <- read.csv("SSI Jan21.csv", header=T)
summary(SSI)
litter.consol$SSI <- SSI$ssi.ba[match(litter.consol$species,SSI$X)]
litter.consol$iv.wet <- SSI$iv.wet[match(litter.consol$species,SSI$X)]
litter.consol$iv.wet[which(is.na(litter.consol$iv.wet))] <- SSI$iv.wet[which(SSI$X=="Madhuca sp.")]
litter.consol$SSI[which(is.na(litter.consol$SSI))] <- SSI$ssi.ba[which(SSI$X=="Madhuca sp.")]

#plot(log(litter.production) ~ iv.wet, data=litter.consol)
#plot(log(litter.production) ~ SSI, data=litter.consol)

litter.consol

########################
# Loglitter linear mod #
########################

# 34 species
unique(litter.consol$species)

summary(litter.consol)
litter.consol$loglitter <- log(litter.consol$litter.production)
litter.consol$logba <- log(litter.consol$ten.ba)
litter.consol$plotwet <- ifelse(litter.consol$plot %in% c("Q3","Q4","Q6","Q9","Q10"), 1, 0)

#############
# CHNS data #
#############

chns <- read.csv("CHNS v3.csv", header=T)
summary(chns)
cn <- with(chns, tapply(Ratio, Species, mean))
litter.consol$CNratio <- cn[match(litter.consol$species, names(cn))]
# number of species represented
sum(ifelse(litter.consol$CNratio != "NA",1,0),na.rm=T)
plot(log(litter.production) ~ CNratio, data=litter.consol, cex=2)

#write.csv(litter.consol, "litter.production w species-plot CAA Feb21.csv")

# Aggregate tree and plot basal area data
treeleaves <- lamina[which(!lamina$species %in% lianas),]
tl.sum <- aggregate(Dry.Mass ~ date1 + date2 + Plot, data=treeleaves, FUN=sum)
plot.type <- unique(lamina[c("Plot","Condition")])
tll <- merge(tl.sum, plot.type)
plot.ba <- with(tree, tapply(ba, plot, sum))
#pt <- c("D","W","D","D","W","W","D","D","W","W")
#boxplot(plot.ba ~ pt)
tll$ba <- plot.ba[match(tll$Plot, names(plot.ba))]

tll$prod <- tll$Dry.Mass / tll$ba*100^2 / 4

summary(tll)
summary(lamina)


se <- function(x) sd(x)/sqrt(length(x))

(all <- with(lamina, tapply(Dry.Mass, Plot, sum)))
(bolian <- with(tll, tapply(Dry.Mass, Plot, sum)))
(bolianba <- with(tll, tapply(prod, Plot, sum)))
# convert to Mg per hectare per year
#all <- 10000*all/4000/1000
#bolian <- 10000*bolian/4000/1000
# convert to g per m^2 per year
all <- all/4
bolian <- bolian/4

pt <- c("D","W","D","D","W","W","D","D","W","W")

#plot(bolian ~ log(plot.ba), pch=16, col=ifelse(pt=="W", col.pal[2], col.pal[1]), cex=3)
#points(all ~ log(plot.ba), pch=21, lwd=2, col=ifelse(pt=="W", col.pal[2], col.pal[1]), cex=3, bg="white")

t.test(all~pt, var.equal=T)
t.test(bolian~pt, var.equal=T)
t.test(bolianba~pt, var.equal=T)

means <- matrix(rep(NA,6), ncol=2, dimnames=list(c("all", "bolian", "bolianba"),c("D","W")))
ses <- matrix(rep(NA,6), ncol=2, dimnames=list(c("all", "bolian", "bolianba"),c("D","W")))

# grand total mean and SE
mean(all)
se(all)

means[1,] <- tapply(all, pt, mean)
means[2,] <- tapply(bolian, pt, mean)
means[3,] <- tapply(bolianba, pt, mean)
ses[1,] <- tapply(all, pt, se)
ses[2,] <- tapply(bolian, pt, se)
ses[3,] <- tapply(bolianba, pt, se)
means
ses
(upp <- means+ses)
(low <- means-ses)

col.pal <- c("#9EC1A3" ,"#40798C",
	"#70A9A1", "#1F363D", "#8D8266")
col.t <- adjustcolor(col.pal, alpha.f=0.7)
#barplot(matrix(rep(1,5)), col=col.t)

#pdf("barplots revised.pdf", height=8, width=14, useDingbats=F)
#jpeg("barplots.jpg", height=8, width=14, units="in", res=300)
par(mfrow=c(1,2), mar=c(4,5.5,2,1.5), mgp=c(3.5,1,0))
bpa <- barplot(means[1:2,], beside=T, ylim=c(0,1000), col=col.pal[c(3,1)],
  names=c("Non-swamp", "Swamp"), cex.lab=1.5, las=1, cex.axis=1.5, cex.names=1.5,
  ylab=expression(paste("Canopy production (g", " year"^"-1", " m"^"-2", ")")))
for(i in 1:4){
  arrows(c(bpa)[i], c(upp[1:2,])[i], c(bpa)[i], c(low[1:2,])[i], angle=90, code=3)}
legend('topright', legend=c("All leaves", "Excluding liana species"), cex=1.5,
  fill=col.pal[c(3,1)])
mtext(side=3, adj=0, line=-1, cex=1.5, text="  a)")
bpb <- barplot(means[3,], beside=T, width=0.5, space=0.75, xlim=c(0,2), ylim=c(0,1000),
  names=c("Non-swamp", "Swamp"), cex.lab=1.5, las=1, cex.axis=1.5, cex.names=1.5, col=col.pal[1],
# ylab=expression(paste("Basal area-normalized canopy production (g", " year"^" -1", " m"^" -4", ")")))
  ylab=expression(paste("Basal area-normalized canopy production (g", " year"^" -1", " m"^" -2", " m"^" -2", ")")))
for(i in 1:2) arrows(bpb[i,1], upp[3,i], bpb[i,1], low[3,i], angle=90, code=3)
mtext(side=3, adj=0, line=-1, cex=1.5, text="  b)")
dev.off()



summary(litter.consol)

TableS1 <- aggregate(cbind(litter.collection, ten.ba) ~ species + code + duration, data=litter.consol, FUN=sum)
TableS1$duration <- round(TableS1$duration / 7, 0)
nplots <- table(litter.consol$species)
TableS1$No.of.plots <- nplots[match(TableS1$species, names(nplots))]
TableS1$litter.production <- round(
	with(TableS1, 
	(litter.collection / (4*No.of.plots)) 	# convert units to g/m2 of plot area
	/ (duration/52))					# per year of survey
					, 1)
TableS1$productivity <- round(TableS1$litter.production / (TableS1$ten.ba/10000), 1)
TableS1$litter.collection <- round(TableS1$litter.collection, 1)
TableS1$ten.ba <- round(TableS1$ten.ba, 1)
TableS1
#write.csv(TableS1, "Table S1 Feb21.csv")






litter.consol <- na.omit(litter.consol)
summary(litter.consol)

library(MuMIn)
summary(w0 <- lm(log(litter.production) ~ 1, weights=logba, data=litter.consol))
summary(w1 <- lm(log(litter.production) ~ SSI, weights=logba, data=litter.consol))
summary(w2 <- lm(log(litter.production) ~ SSI + I(SSI^2), weights=logba, data=litter.consol))
summary(w1s <- lm(log(litter.production) ~ SSI.stem, weights=logba, data=litter.consol))
summary(w2s <- lm(log(litter.production) ~ SSI.stem + I(SSI.stem^2), weights=logba, data=litter.consol))

summary(m0 <- lm(log(litter.production) ~ 1, data=litter.consol))
summary(m1 <- lm(log(litter.production) ~ SSI, data=litter.consol))
summary(m2 <- lm(log(litter.production) ~ SSI + I(SSI^2), data=litter.consol))
summary(m1s <- lm(log(litter.production) ~ SSI.stem, data=litter.consol))
summary(m2s <- lm(log(litter.production) ~ SSI.stem + I(SSI.stem^2), data=litter.consol))

summary(b0 <- lm(log(litter.production) ~ logba, data=litter.consol))
summary(b1 <- lm(log(litter.production) ~ logba + SSI, data=litter.consol))
summary(b2 <- lm(log(litter.production) ~ logba + SSI + I(SSI^2), data=litter.consol))
summary(b3 <- lm(log(litter.production) ~ logba*SSI, data=litter.consol))
summary(b4 <- lm(log(litter.production) ~ logba*SSI + I(SSI^2), data=litter.consol))
summary(b5 <- lm(log(litter.production) ~ logba*SSI + logba*I(SSI^2), data=litter.consol))
summary(b1s <- lm(log(litter.production) ~ logba + SSI.stem, data=litter.consol))
summary(b2s <- lm(log(litter.production) ~ logba + SSI.stem + I(SSI.stem^2), data=litter.consol))
summary(b3s <- lm(log(litter.production) ~ logba*SSI.stem, data=litter.consol))
summary(b4s <- lm(log(litter.production) ~ logba*SSI.stem + I(SSI.stem^2), data=litter.consol))
summary(b5s <- lm(log(litter.production) ~ logba*SSI.stem + logba*I(SSI.stem^2), data=litter.consol))

(mods <- model.sel(w0,w1,w2,m0,m1,m2,b0,b1,b2,b3,b4,b5,
	w1s,w2s,m1s,m2s,b1s,b2s,b3s,b4s,b5s))

x <- seq(-1,1,0.01)
pred <- predict(w1, newdata=data.frame(SSI=x), se.fit=T)



par(mar=c(6,6,1.5,1.5),mgp=c(4,1,0))
plot(log(litter.production) ~ SSI, type="n", xlim=c(-1.2,1.2), cex.lab=1.5, cex.axis=1.5,
	xlab="Swamp specialization index", data=litter.consol, yaxt="n",
	ylab=expression(paste("Species litter production constant (g", " m"^" -2", " year"^" -1", "  m"^" -2", ")")))
axis(side=2, at=log(c(10,50,100,500,1000,5000)), labels=c(10,50,100,500,1000,5000), las=1, cex.axis=1.5)
polygon(c(pred$fit+1.96*pred$se.fit, rev(pred$fit-1.96*pred$se.fit)) ~ c(x, rev(x)), border=F, col="grey")
lines(pred$fit ~ x, lwd=5, col="white")
text(log(litter.production) ~ SSI, labels=species, cex=(logba-3)/3, font=3, data=litter.consol)
axis(side=1, at=c(-1,1), mgp=c(4,2.5,0), labels=c("Dry specialist","Wet specialist"), cex.axis=1.5)



