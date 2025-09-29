library(tidyverse)

PJtraits <- read.csv("./data/PJ leaf dat.csv", header=T)

HJtraits <- read.csv("../ecophysio_traits/raw_data/leaf soft traits.csv", header=T) %>%
  mutate(ldmc = as.numeric(ldmc),
         SLA = as.numeric(SLA),
         thickness = as.numeric(thickness),
         species = toupper(substr(indiv, 0, 3)))

traits <- rbind(
  HJtraits %>%
    group_by(species) %>%
    summarise(sla = mean(SLA, na.rm = TRUE),
              lt = mean(thickness, na.rm = TRUE),
              ldmc = mean(ldmc, na.rm = TRUE)),
  PJtraits %>%
    mutate(meanlt = mean(c(MESOPHYLL_THICKNESS_1_mm,
                           MESOPHYLL_THICKNESS_2_mm,
                           MESOPHYLL_THICKNESS_3_mm))) %>%
    rename(species = SPECIES) %>%
    group_by(species) %>%
    summarise(sla = mean(SLA_cm2.g.1, na.rm = TRUE),
              lt = mean(meanlt, na.rm = TRUE),
              ldmc = mean(LDMC, na.rm = TRUE))
)


PJsla <- with(PJtraits, tapply(SLA_cm2.g.1, SPECIES, mean, na.rm=T))
PJtraits$meanlt <- apply(PJtraits[,6:8],1,mean,na.rm=T)
PJlt <- with(PJtraits, tapply(meanlt, SPECIES, mean))
PJldmc <- with(PJtraits, tapply(LDMC, SPECIES, mean, na.rm=T))


# number of individuals represented in the samples
(PJsample <- table(PJtraits$SPECIES)/5)
(HJsample <- table(HJtraits$species)/5)

traitsonly <- data.frame(code=c(names(HJsla),names(PJsla)), SLA=c(HJsla,PJsla), LDMC=c(HJldmc, PJldmc), LT=c(HJlt,PJlt))
cnr <- with(litter, tapply(CNratio, code, mean))
traitsonly <- merge(traitsonly, data.frame(code=names(cnr),CNratio=cnr), by="code")
leaf.PCA <- prcomp(traitsonly[,2:5], scale=T)
summary(leaf.PCA)
rownames(leaf.PCA$rotation)[4] <- "C:N"
biplot(leaf.PCA)

###################################
# Merge trait and production data #
###################################

litter <- read.csv("./data/litter.production w species-plot CAA Feb21.csv", header=T) %>%
  select(-X)

litter %>%
  group_by(code) %>%
  summarise(CNratio = mean(CNratio)) %>%
  rename(species = code) %>%
  full_join(traits, by = "species") %>%
  View()

litter <- rbind(
  merge(x=litter, y=data.frame(code=names(HJsla),SLA=HJsla), by="code"),
  merge(x=litter, y=data.frame(code=names(PJsla),SLA=PJsla), by="code") )
litter <- rbind(
  merge(x=litter, y=data.frame(code=names(HJlt),LT=HJlt), by="code"),
  merge(x=litter, y=data.frame(code=names(PJlt),LT=PJlt), by="code") )
litter <- rbind(
  merge(x=litter, y=data.frame(code=names(HJldmc),LDMC=HJldmc), by="code"),
  merge(x=litter, y=data.frame(code=names(PJldmc),LDMC=PJldmc), by="code") )
pcax <- data.frame(leaf.PCA$x)
pcax$code <- traitsonly$code
litter <- merge(x=litter, y=pcax, by="code")

summary(litter)

nrow(litter)
unique(litter$species)
# assign CTE the mean value for now
litter$CNratio[which(litter$code=="CTE")] <- mean(traitsonly$CNratio, na.rm=T)
hist(litter$loglitter <- log(litter$litter.production))
hist(litter$logba <- log(litter$ten.ba))

# weight individual points by the inverse of the number of "pseudo replicates" for that species
(sp.weights <- 1/table(litter$species))
litter$weights <- as.numeric(sp.weights[match(litter$species, names(sp.weights))])
with(litter, tapply(weights, species, sum))

# include plot type variable
litter$plotwet <- ifelse(litter$plot %in% c("Q3","Q4","Q6","Q9","Q10"), 1, 0)
boxplot(SSI~plotwet, data=litter)

#rsi <- read.csv("RSI 2018.csv")
#litter$SSI <- rsi$rsi.ba[match(litter$species, rsi$X)]
#litter$SSI.stem <- rsi$rsi.stem[match(litter$species, rsi$X)]

####################
# COMPONENT MODELS #
####################

library(piecewiseSEM)
library(lme4)

summary(litter)
litter$SLA <- litter$SLA/100

spweights <- table(litter$species)^-1
litter$weights <- spweights[match(litter$species, names(spweights))]

pairs.cor(litter[c("logba","loglitter","SSI","SLA","LT","LDMC","CNratio","PC1","PC2","plotwet")])

# iv.dry not very useful, iv.wet must transform
litter$iv.wet <- sqrt(litter$iv.wet)
summary(litter)

# Leaf traits
summary(S10 <- lm(SLA ~ 1, weights=weights, data=litter))
summary(S1a <- lm(SLA ~ SSI, weights=weights, data=litter))
summary(S1b <- lm(SLA ~ SSI.stem, weights=weights, data=litter))
summary(S1c <- lm(SLA ~ iv.wet, weights=weights, data=litter))

summary(S20 <- lm(LT ~ 1, weights=weights, data=litter))
summary(S2a <- lm(LT ~ SSI, weights=weights, data=litter))
summary(S2b <- lm(LT ~ SSI.stem, weights=weights, data=litter))
summary(S2c <- lm(LT ~ iv.wet, weights=weights, data=litter))

summary(S30 <- lm(LDMC ~ 1, weights=weights, data=litter))
summary(S3a <- lm(LDMC ~ SSI, weights=weights, data=litter))
summary(S3b <- lm(LDMC ~ SSI.stem, data=litter))
summary(S3c <- lm(LDMC ~ iv.wet, weights=weights, data=litter))

summary(S40 <- lm(CNratio ~ 1, weights=weights, data=litter))
summary(S4a <- lm(CNratio ~ SSI, weights=weights, data=litter))
summary(S4b <- lm(CNratio ~ SSI.stem, weights=weights, data=litter))
summary(S4c <- lm(CNratio ~ iv.wet, weights=weights, data=litter))

summary(S50 <- lm(PC1 ~ 1, weights=weights, data=litter))
summary(S5a <- lm(PC1 ~ SSI, weights=weights, data=litter))
summary(S5b <- lm(PC1 ~ SSI.stem, weights=weights, data=litter))
summary(S5c <- lm(PC1 ~ iv.wet, weights=weights, data=litter))

summary(S60 <- lm(PC2 ~ 1, weights=weights, data=litter))
summary(S6a <- lm(PC2 ~ SSI, weights=weights, data=litter))
summary(S6b <- lm(PC2 ~ SSI.stem, weights=weights, data=litter))
summary(S6c <- lm(PC2 ~ iv.wet, weights=weights, data=litter))

# litter production
summary(L0 <- lmer(loglitter ~ logba + (1|plot) + (logba|species), data=litter))

summary(L1 <- lmer(loglitter ~ logba + SLA + (1|plot) + (logba|species), data=litter))
summary(L2 <- lmer(loglitter ~ logba + LT + (1|plot) + (logba|species), data=litter))
summary(L3 <- lmer(loglitter ~ logba + LDMC + (1|plot) + (logba|species), data=litter))
summary(L4 <- lmer(loglitter ~ logba + CNratio + (1|plot) + (logba|species), data=litter))
summary(LP <- lmer(loglitter ~ logba + PC1 + (1|plot) + (logba|species), data=litter))
summary(LQ <- lmer(loglitter ~ logba + PC2 + (1|plot) + (logba|species), data=litter))

summary(L5 <- lmer(loglitter ~ logba + SSI + (1|plot) + (logba|species), data=litter))
summary(L6 <- lmer(loglitter ~ logba + SSI.stem + (1|plot) + (logba|species), data=litter))
summary(L7 <- lmer(loglitter ~ logba + iv.wet + (1|plot) + (logba|species), data=litter))

summary(L8.1a <- lmer(loglitter ~ logba + SSI + SLA + (1|plot) + (logba|species), data=litter))
summary(L8.1b <- lmer(loglitter ~ logba + SSI.stem + SLA + (1|plot) + (logba|species), data=litter))
summary(L8.1c <- lmer(loglitter ~ logba + iv.wet + SLA + (1|plot) + (logba|species), data=litter))

summary(L8.2a <- lmer(loglitter ~ logba + SSI + LT + (1|plot) + (logba|species), data=litter))
summary(L8.2b <- lmer(loglitter ~ logba + SSI.stem + LT + (1|plot) + (logba|species), data=litter))
summary(L8.2c <- lmer(loglitter ~ logba + iv.wet + LT + (1|plot) + (logba|species), data=litter))

summary(L8.3a <- lmer(loglitter ~ logba + SSI + LDMC + (1|plot) + (logba|species), data=litter))
summary(L8.3b <- lmer(loglitter ~ logba + SSI.stem + LDMC + (1|plot) + (logba|species), data=litter))
summary(L8.3c <- lmer(loglitter ~ logba + iv.wet + LDMC + (1|plot) + (logba|species), data=litter))

summary(L8.4a <- lmer(loglitter ~ logba + SSI + CNratio + (1|plot) + (logba|species), data=litter))
summary(L8.4b <- lmer(loglitter ~ logba + SSI.stem + CNratio + (1|plot) + (logba|species), data=litter))
summary(L8.4c <- lmer(loglitter ~ logba + iv.wet + CNratio + (1|plot) + (logba|species), data=litter))

summary(L8.5a <- lmer(loglitter ~ logba + SSI + PC1 + (1|plot) + (logba|species), data=litter))
summary(L8.5b <- lmer(loglitter ~ logba + SSI.stem + PC1 + (1|plot) + (logba|species), data=litter))
summary(L8.5c <- lmer(loglitter ~ logba + iv.wet + PC1 + (1|plot) + (logba|species), data=litter))

summary(L8.6a <- lmer(loglitter ~ logba + SSI + PC2 + (1|plot) + (logba|species), data=litter))
summary(L8.6b <- lmer(loglitter ~ logba + SSI.stem + PC2 + (1|plot) + (logba|species), data=litter))
summary(L8.6c <- lmer(loglitter ~ logba + iv.wet + PC2 + (1|plot) + (logba|species), data=litter))

# with plot type
summary(L0t <- lmer(loglitter ~ logba + (1|plot) + plotwet + (logba|species), data=litter))

summary(L1t <- lmer(loglitter ~ logba + SLA + plotwet + (1|plot) + (logba|species), data=litter))
summary(L2t <- lmer(loglitter ~ logba + LT + plotwet + (1|plot) + (logba|species), data=litter))
summary(L3t <- lmer(loglitter ~ logba + LDMC + plotwet + (1|plot) + (logba|species), data=litter))
summary(L4t <- lmer(loglitter ~ logba + CNratio + plotwet + (1|plot) + (logba|species), data=litter))
summary(LPt <- lmer(loglitter ~ logba + PC1 + plotwet + (1|plot) + (logba|species), data=litter))
summary(LQt <- lmer(loglitter ~ logba + PC2 + plotwet + (1|plot) + (logba|species), data=litter))

summary(L5t <- lmer(loglitter ~ logba + SSI + plotwet + (1|plot) + (logba|species), data=litter))
summary(L6t <- lmer(loglitter ~ logba + SSI.stem + plotwet + (1|plot) + (logba|species), data=litter))
summary(L7t <- lmer(loglitter ~ logba + iv.wet + plotwet + (1|plot) + (logba|species), data=litter))

summary(L8.1at <- lmer(loglitter ~ logba + SSI + SLA + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.1bt <- lmer(loglitter ~ logba + SSI.stem + SLA + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.1ct <- lmer(loglitter ~ logba + iv.wet + SLA + plotwet + (1|plot) + (logba|species), data=litter))

summary(L8.2at <- lmer(loglitter ~ logba + SSI + LT + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.2bt <- lmer(loglitter ~ logba + SSI.stem + LT + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.2ct <- lmer(loglitter ~ logba + iv.wet + LT + plotwet + (1|plot) + (logba|species), data=litter))

summary(L8.3at <- lmer(loglitter ~ logba + SSI + LDMC + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.3bt <- lmer(loglitter ~ logba + SSI.stem + LDMC + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.3ct <- lmer(loglitter ~ logba + iv.wet + LDMC + plotwet + (1|plot) + (logba|species), data=litter))

summary(L8.4at <- lmer(loglitter ~ logba + SSI + CNratio + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.4bt <- lmer(loglitter ~ logba + SSI.stem + CNratio + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.4ct <- lmer(loglitter ~ logba + iv.wet + CNratio + plotwet + (1|plot) + (logba|species), data=litter))

summary(L8.5at <- lmer(loglitter ~ logba + SSI + PC1 + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.5bt <- lmer(loglitter ~ logba + SSI.stem + PC1 + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.5ct <- lmer(loglitter ~ logba + iv.wet + PC1 + plotwet + (1|plot) + (logba|species), data=litter))

summary(L8.6at <- lmer(loglitter ~ logba + SSI + PC2 + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.6bt <- lmer(loglitter ~ logba + SSI.stem + PC2 + plotwet + (1|plot) + (logba|species), data=litter))
summary(L8.6ct <- lmer(loglitter ~ logba + iv.wet + PC2 + plotwet + (1|plot) + (logba|species), data=litter))

# distribution
summary(D0a <- lm(SSI ~ 1, weights=weights, data=litter))
summary(D0b <- lm(SSI.stem ~ 1, weights=weights, data=litter))
summary(D0c <- lm(iv.wet ~ 1, weights=weights, data=litter))

summary(D1a <- lm(SSI ~ SLA, weights=weights, data=litter))
summary(D1b <- lm(SSI.stem ~ SLA, weights=weights, data=litter))
summary(D1c <- lm(iv.wet ~ SLA, weights=weights, data=litter))
summary(D2a <- lm(SSI ~ LT, weights=weights, data=litter))
summary(D2b <- lm(SSI.stem ~ LT, weights=weights, data=litter))
summary(D2c <- lm(iv.wet ~ LT, weights=weights, data=litter))
summary(D3a <- lm(SSI ~ LDMC, weights=weights, data=litter))
summary(D3b <- lm(SSI.stem ~ LDMC, weights=weights, data=litter))
summary(D3c <- lm(iv.wet ~ LDMC, weights=weights, data=litter))
summary(D4a <- lm(SSI ~ CNratio, weights=weights, data=litter))
summary(D4b <- lm(SSI.stem ~ CNratio, weights=weights, data=litter))
summary(D4c <- lm(iv.wet ~ CNratio, weights=weights, data=litter))
summary(D5a <- lm(SSI ~ PC1, weights=weights, data=litter))
summary(D5b <- lm(SSI.stem ~ PC1, weights=weights, data=litter))
summary(D5c <- lm(iv.wet ~ PC1, weights=weights, data=litter))
summary(D6a <- lm(SSI ~ PC2, weights=weights, data=litter))
summary(D6b <- lm(SSI.stem ~ PC2, weights=weights, data=litter))
summary(D6c <- lm(iv.wet ~ PC2, weights=weights, data=litter))

###################
# PATH MODEL LIST #
###################

mods <- list()

### triangular models
## A: leaf traits determines both HS and CP, HS has an additional effect on CP
# SLA
mods[[1]] <- psem(L8.1c, D1c)
# LT
mods[[length(mods)+1]] <- psem(L8.2c, D2c)
# LDMC
mods[[length(mods)+1]] <- psem(L8.3c, D3c)
# CN ratio
mods[[length(mods)+1]] <- psem(L8.4c, D4c)
# PC1
mods[[length(mods)+1]] <- psem(L8.5c, D5c)
# PC2
mods[[length(mods)+1]] <- psem(L8.6c, D6c)

## B: HS determines leaf traits and CP, leaf traits have an additional effect on CP
# SLA
mods[[length(mods)+1]] <- psem(L8.1c, S1c)
# LT
mods[[length(mods)+1]] <- psem(L8.2c, S2c)
# LDMC
mods[[length(mods)+1]] <- psem(L8.3c, S3c)
# CN ratio
mods[[length(mods)+1]] <- psem(L8.4c, S4c)
# PC1
mods[[length(mods)+1]] <- psem(L8.5c, S5c)
# PC2
mods[[length(mods)+1]] <- psem(L8.6c, S6c)

# Biological hypothesis code for triangular models
BH <- c(rep("i", 6), rep("ii", 6))

### chain models
## A: leaf traits determine specialization, which alters litter productivity
# SLA
# i~13th
mods[[length(mods)+1]] <- psem(D1c, L7)
# LT
mods[[length(mods)+1]] <- psem(D2c, L7)
# LDMC
mods[[length(mods)+1]] <- psem(D3c, L7)
# CN ratio
mods[[length(mods)+1]] <- psem(D4c, L7)
# PC1
mods[[length(mods)+1]] <- psem(D5c, L7)
# PC2
mods[[length(mods)+1]] <- psem(D6c, L7)

## B: specialization constrains leaf traits, which alter litter productivity
# i~19th
mods[[length(mods)+1]] <- psem(S1c, L1)
mods[[length(mods)+1]] <- psem(S2c, L2)
mods[[length(mods)+1]] <- psem(S3c, L3)
mods[[length(mods)+1]] <- psem(S4c, L4)
mods[[length(mods)+1]] <- psem(S5c, LP)
mods[[length(mods)+1]] <- psem(S6c, LQ)

# Biological hypothesis code for chain models
BH <- c(BH, rep("vi", 6), rep("vii", 6))

### LAMBDA model
# SLA
mods[[length(mods)+1]] <- psem(L8.1c, S10, D0c)
# LT
mods[[length(mods)+1]] <- psem(L8.2c, S20, D0c)
# LDMC
mods[[length(mods)+1]] <- psem(L8.3c, S30, D0c)
# CN ratio
mods[[length(mods)+1]] <- psem(L8.4c, S40, D0c)
# PC1
mods[[length(mods)+1]] <- psem(L8.5c, S50, D0c)
# PC2
mods[[length(mods)+1]] <- psem(L8.6c, S60, D0c)

# Biological hypothesis code for lambda model
BH <- c(BH, rep("iii", 6))

### Single cause (null) models
## A: specialization determines litter production
mods[[length(mods)+1]] <- psem(L7, D0c, S10)
mods[[length(mods)+1]] <- psem(L7, D0c, S20)
mods[[length(mods)+1]] <- psem(L7, D0c, S30)
mods[[length(mods)+1]] <- psem(L7, D0c, S40)
mods[[length(mods)+1]] <- psem(L7, D0c, S50)
mods[[length(mods)+1]] <- psem(L7, D0c, S60)

## B: leaf traits determine distribution
mods[[length(mods)+1]] <- psem(D1c, S10, L0)
mods[[length(mods)+1]] <- psem(D2c, S20, L0)
mods[[length(mods)+1]] <- psem(D3c, S30, L0)
mods[[length(mods)+1]] <- psem(D4c, S40, L0)
mods[[length(mods)+1]] <- psem(D5c, S50, L0)
mods[[length(mods)+1]] <- psem(D6c, S60, L0)

## C: habitat specialization constrains leaf traits
mods[[length(mods)+1]] <- psem(S1c, D0c, L0)
mods[[length(mods)+1]] <- psem(S2c, D0c, L0)
mods[[length(mods)+1]] <- psem(S3c, D0c, L0)
mods[[length(mods)+1]] <- psem(S4c, D0c, L0)
mods[[length(mods)+1]] <- psem(S5c, D0c, L0)
mods[[length(mods)+1]] <- psem(S6c, D0c, L0)

## D: leaf traits determine productivity
mods[[length(mods)+1]] <- psem(L1, S10, D0c)
mods[[length(mods)+1]] <- psem(L2, S20, D0c)
mods[[length(mods)+1]] <- psem(L3, S30, D0c)
mods[[length(mods)+1]] <- psem(L4, S40, D0c)
mods[[length(mods)+1]] <- psem(LP, S50, D0c)
mods[[length(mods)+1]] <- psem(LQ, S60, D0c)

# Biological hypothesis code for single cause (null) models
BH <- c(BH, rep("iv",6), rep("viii",6), rep("ix",6), rep("v",6))
#BH <- c(BH, rep("viii",6), rep("ix",6))
length(BH); length(mods) #42

### plus plot type models (exact repeat of all the previous models)

### triangular models
## A: leaf traits determines both HS and CP, HS has an additional effect on CP
mods[[length(mods)+1]] <- psem(L8.1ct, S1c)
mods[[length(mods)+1]] <- psem(L8.2ct, D2c)
mods[[length(mods)+1]] <- psem(L8.3ct, D3c)
mods[[length(mods)+1]] <- psem(L8.4ct, D4c)
mods[[length(mods)+1]] <- psem(L8.5ct, D5c)
mods[[length(mods)+1]] <- psem(L8.6ct, D6c)

## B: HS determines leaf traits and CP, leaf traits have an additional effect on CP
mods[[length(mods)+1]] <- psem(L8.1ct, S1c)
mods[[length(mods)+1]] <- psem(L8.2ct, S2c)
mods[[length(mods)+1]] <- psem(L8.3ct, S3c)
mods[[length(mods)+1]] <- psem(L8.4ct, S4c)
mods[[length(mods)+1]] <- psem(L8.5ct, S5c)
mods[[length(mods)+1]] <- psem(L8.6ct, S6c)

# Biological hypothesis code for triangular models
BH <- c(BH, rep("i", 6), rep("ii", 6))

### chain models
## A: leaf traits determine specialization, which alters litter productivity
# i~55th
mods[[length(mods)+1]] <- psem(D1c, L7t)
mods[[length(mods)+1]] <- psem(D2c, L7t)
mods[[length(mods)+1]] <- psem(D3c, L7t)
mods[[length(mods)+1]] <- psem(D4c, L7t)
mods[[length(mods)+1]] <- psem(D5c, L7t)
mods[[length(mods)+1]] <- psem(D6c, L7t)

## B: specialization constrains leaf traits, which alter litter productivity
# i~61st
mods[[length(mods)+1]] <- psem(S1c, L1t)
mods[[length(mods)+1]] <- psem(S2c, L2t)
mods[[length(mods)+1]] <- psem(S3c, L3t)
mods[[length(mods)+1]] <- psem(S4c, L4t)
mods[[length(mods)+1]] <- psem(S5c, LPt)
mods[[length(mods)+1]] <- psem(S6c, LQt)

# Biological hypothesis code for chain models
BH <- c(BH, rep("vi", 6), rep("vii", 6))

### LAMBDA model
mods[[length(mods)+1]] <- psem(L8.1ct, S10, D0c)
mods[[length(mods)+1]] <- psem(L8.2ct, S20, D0c)
mods[[length(mods)+1]] <- psem(L8.3ct, S30, D0c)
mods[[length(mods)+1]] <- psem(L8.4ct, S40, D0c)
mods[[length(mods)+1]] <- psem(L8.5ct, S50, D0c)
mods[[length(mods)+1]] <- psem(L8.6ct, S60, D0c)

# Biological hypothesis code for lambda model
BH <- c(BH, rep("iii", 6))

### Single cause (null) models
## A: specialization determines litter production
mods[[length(mods)+1]] <- psem(L7t, D0c, S10)
mods[[length(mods)+1]] <- psem(L7t, D0c, S20)
mods[[length(mods)+1]] <- psem(L7t, D0c, S30)
mods[[length(mods)+1]] <- psem(L7t, D0c, S40)
mods[[length(mods)+1]] <- psem(L7t, D0c, S50)
mods[[length(mods)+1]] <- psem(L7t, D0c, S60)

## B: leaf traits determine distribution
mods[[length(mods)+1]] <- psem(D1c, S10, L0t)
mods[[length(mods)+1]] <- psem(D2c, S20, L0t)
mods[[length(mods)+1]] <- psem(D3c, S30, L0t)
mods[[length(mods)+1]] <- psem(D4c, S40, L0t)
mods[[length(mods)+1]] <- psem(D5c, S50, L0t)
mods[[length(mods)+1]] <- psem(D6c, S60, L0t)

## C: habitat specialization constrains leaf traits
mods[[length(mods)+1]] <- psem(S1c, D0c, L0t)
mods[[length(mods)+1]] <- psem(S2c, D0c, L0t)
mods[[length(mods)+1]] <- psem(S3c, D0c, L0t)
mods[[length(mods)+1]] <- psem(S4c, D0c, L0t)
mods[[length(mods)+1]] <- psem(S5c, D0c, L0t)
mods[[length(mods)+1]] <- psem(S6c, D0c, L0t)

## D: leaf traits determine productivity
mods[[length(mods)+1]] <- psem(L1t, S10, D0c)
mods[[length(mods)+1]] <- psem(L2t, S20, D0c)
mods[[length(mods)+1]] <- psem(L3t, S30, D0c)
mods[[length(mods)+1]] <- psem(L4t, S40, D0c)
mods[[length(mods)+1]] <- psem(LPt, S50, D0c)
mods[[length(mods)+1]] <- psem(LQt, S60, D0c)

# Biological hypothesis code for single cause (null) models
BH <- c(BH, rep("iv",6), rep("viii",6), rep("ix",6), rep("v",6))
#BH <- c(BH, rep("viii",6), rep("ix",6))

length(BH); length(mods) #108

###################
# MODEL SELECTION #
###################

AIC.list<-rep(NA,length(mods))
for(i in 1:length(mods)){ tryCatch({	
	AIC.list[i] <- AIC(mods[[i]],aicc=T)
	},error=function(e){})}
min(AIC.list, na.rm=T)
dAIC <- AIC.list-min(AIC.list, na.rm=T)
# rank
#(best.index <- order(dAIC)[1:30])
(ranked.index <- order(dAIC))
#(dAIC.ranked <- round(dAIC[best.index],2))
(dAIC.ranked <- round(dAIC[ranked.index],2))

summary(mods[[ ranked.index[1] ]], .progressBar=F) # SLA: -0.2905*-1.6676 = 0.4844
summary(mods[[ ranked.index[2] ]], .progressBar=F) # PC1: 0.0707*-1.6676 = -0.1179
summary(mods[[ ranked.index[3] ]], .progressBar=F) # BH i. SLA direct eff = 0.1496
summary(mods[[ ranked.index[4] ]], .progressBar=F) # BH iv. LDMC direct eff = 2.2548
summary(mods[[ ranked.index[5] ]], .progressBar=F) # LT: 1.1955*-1.6676 = -1.9936
summary(mods[[ ranked.index[6] ]], .progressBar=F) # BH i. PC1 direct eff = 0.1336
summary(mods[[ ranked.index[7] ]], .progressBar=F) 
summary(mods[[ ranked.index[8] ]], .progressBar=F) # CNratio: 0.0088*-1.6676 = -0.0147
summary(mods[[ ranked.index[9] ]], .progressBar=F) 
summary(mods[[ ranked.index[10] ]], .progressBar=F)
summary(mods[[ ranked.index[11] ]], .progressBar=F)
summary(mods[[ ranked.index[12] ]], .progressBar=F)
summary(mods[[ ranked.index[13] ]], .progressBar=F)
summary(mods[[ ranked.index[14] ]], .progressBar=F)
summary(mods[[ ranked.index[15] ]], .progressBar=F)
summary(mods[[ ranked.index[16] ]], .progressBar=F)
summary(mods[[ ranked.index[17] ]], .progressBar=F)
summary(mods[[ ranked.index[25] ]], .progressBar=F)
rsquared(mods[[ ranked.index[2] ]])

### SUM OF WEIGHTS OF BIOLOGICAL HYPOTHESES ###
unranked.AICC <- round(AIC.list,2)
unranked.RL <- round(exp(-0.5*dAIC),5)
unranked.weights <- round(unranked.RL/sum(unranked.RL),3)
SoW <- tapply(unranked.weights, BH, sum)
SoW

### AIC TABLE ###

AICC <- round(AIC.list[ranked.index],2)
RL <- round(exp(-0.5*dAIC.ranked),5)
weights <- round(RL/sum(RL, na.rm=T),3)
BH.ranked <- BH[ranked.index]

final.table <- data.frame(
	formula = NA, 
	Biol.Hypothesis = BH.ranked,
	R2 = NA,
	Fisher.C = NA,
	AICc = AICC, dAICc = dAIC.ranked, weights = weights)

for(i in 1:nrow(final.table)){
	final.table$R2[i] <-
	paste(round(100*rsquared(mods[[ ranked.index[i] ]])$Marginal,2), collapse="\n")
	}
for(i in 1:nrow(final.table)){
	final.table$formula[i] <-
	summary(mods[[ ranked.index[i] ]])$call
	}
for(i in 1:nrow(final.table)){
	final.table$Fisher.C[i] <-
	as.numeric(summary(mods[[ ranked.index[i] ]])$Cstat[1])
	}
(final.table <- na.omit(final.table))
#write.csv(final.table, "table (iv.wet models) Feb21.csv")

library(MuMIn)
mod2 <- lmer(loglitter ~ logba + iv.wet + plotwet + (1|plot) + (1|species), data=litter)
mod1b <- lmer(loglitter ~ logba + plotwet + (1|plot) + (1|species), data=litter)
mod1a <- lmer(loglitter ~ logba + iv.wet + (1|plot) + (1|species), data=litter)
mod0 <- lmer(loglitter ~ logba + (1|plot) + (1|species), data=litter)
AIC(mod2, mod1b, mod1a, mod0)
AICc(mod2, mod1b, mod1a, mod0)

summary(mod1a)

plot(loglitter ~ iv.wet, data=litter)
plot(resid(mod0) ~ litter$iv.wet)

plot(iv.wet ~ LT, data=litter)





#########
# Plots #
#########

sis <- data.frame(unique(cbind(litter$SLA, litter$iv.wet, litter$code)))
names(sis) <- c("SLA", "iv.wet", "code")
sis$SLA <- as.numeric(sis$SLA)*100

#pdf("Fig 3 new.pdf", width=8, height=15, useDingbats=F)
#jpeg("Fig 3 new.jpg", width=8, height=15, res=300, units="in")
par(mfrow=c(2,1), mar=c(5,5.5,2,2), mgp=c(3.5,1,0))
plot(leaf.PCA$x[,2] ~ leaf.PCA$x[,1], type="n",
	ylab="PC 2 (38.7% of variance)", xlab="PC 1 (48.3% of variance)",
	cex.lab=2, cex.axis=1.5, las=1, xlim=c(-4, 4), ylim=c(-3.5, 2.6))
abline(v=0, lty=2, lwd=2)
abline(h=0, lty=2, lwd=2)
text(leaf.PCA$x[,1], jitter(leaf.PCA$x[,2],350), labels=traitsonly$code, cex=1.5)
for(i in 1:4){
arrows(0, 0, 3*leaf.PCA$rotation[i,1], 3*leaf.PCA$rotation[i,2], col="forestgreen", lwd=3)}
text(4*leaf.PCA$rotation[,1], 3.5*leaf.PCA$rotation[,2], col="forestgreen", cex=2,
	labels=rownames(leaf.PCA$rotation))
plot(iv.wet ~ SLA, data=sis, type="n", xlim=c(60,230),
	xlab=expression(paste("SLA (","cm"^"2", "g"^"-1",")")), ylab="Habitat specialization (swamp indicator value)",
	cex.lab=2, cex.axis=1.5, las=1)
text(sis$SLA, as.numeric(sis$iv.wet), labels=sis$code, cex=1.5)
dev.off()






### PSEM MODEL LIST FOR SSI

mods <- list()

### triangular models
## A: leaf traits determines both HS and CP, HS has an additional effect on CP
# SLA
mods[[1]] <- psem(L8.1a, D1a)
# LT
mods[[length(mods)+1]] <- psem(L8.2a, D2a)
# LDMC
mods[[length(mods)+1]] <- psem(L8.3a, D3a)
# CN ratio
mods[[length(mods)+1]] <- psem(L8.4a, D4a)
# PC1
mods[[length(mods)+1]] <- psem(L8.5a, D5a)
# PC2
mods[[length(mods)+1]] <- psem(L8.6a, D6a)

## B: HS determines leaf traits and CP, leaf traits have an additional effect on CP
# SLA
mods[[length(mods)+1]] <- psem(L8.1a, S1a)
# LT
mods[[length(mods)+1]] <- psem(L8.2a, S2a)
# LDMC
mods[[length(mods)+1]] <- psem(L8.3a, S3a)
# CN ratio
mods[[length(mods)+1]] <- psem(L8.4a, S4a)
# PC1
mods[[length(mods)+1]] <- psem(L8.5a, S5a)
# PC2
mods[[length(mods)+1]] <- psem(L8.6a, S6a)

# Biological hypothesis code for triangular models
BH <- c(rep("i", 6), rep("ii", 6))

### chain models
## A: leaf traits determine specialization, which alters litter productivity
# SLA
# i~13th
mods[[length(mods)+1]] <- psem(D1a, L5)
# LT
mods[[length(mods)+1]] <- psem(D2a, L5)
# LDMC
mods[[length(mods)+1]] <- psem(D3a, L5)
# CN ratio
mods[[length(mods)+1]] <- psem(D4a, L5)
# PC1
mods[[length(mods)+1]] <- psem(D5a, L5)
# PC2
mods[[length(mods)+1]] <- psem(D6a, L5)

## B: specialization constrains leaf traits, which alter litter productivity
# i~19th
mods[[length(mods)+1]] <- psem(S1a, L1)
mods[[length(mods)+1]] <- psem(S2a, L2)
mods[[length(mods)+1]] <- psem(S3a, L3)
mods[[length(mods)+1]] <- psem(S4a, L4)
mods[[length(mods)+1]] <- psem(S5a, LP)
mods[[length(mods)+1]] <- psem(S6a, LQ)

# Biological hypothesis code for chain models
BH <- c(BH, rep("vi", 6), rep("vii", 6))

### LAMBDA model
# SLA
mods[[length(mods)+1]] <- psem(L8.1a, S10, D0a)
# LT
mods[[length(mods)+1]] <- psem(L8.2a, S20, D0a)
# LDMC
mods[[length(mods)+1]] <- psem(L8.3a, S30, D0a)
# CN ratio
mods[[length(mods)+1]] <- psem(L8.4a, S40, D0a)
# PC1
mods[[length(mods)+1]] <- psem(L8.5a, S50, D0a)
# PC2
mods[[length(mods)+1]] <- psem(L8.6a, S60, D0a)

# Biological hypothesis code for lambda model
BH <- c(BH, rep("iii", 6))

### Single cause (null) models
## A: specialization determines litter production
mods[[length(mods)+1]] <- psem(L5, D0a, S10)
mods[[length(mods)+1]] <- psem(L5, D0a, S20)
mods[[length(mods)+1]] <- psem(L5, D0a, S30)
mods[[length(mods)+1]] <- psem(L5, D0a, S40)
mods[[length(mods)+1]] <- psem(L5, D0a, S50)
mods[[length(mods)+1]] <- psem(L5, D0a, S60)

## B: leaf traits determine distribution
mods[[length(mods)+1]] <- psem(D1a, S10, L0)
mods[[length(mods)+1]] <- psem(D2a, S20, L0)
mods[[length(mods)+1]] <- psem(D3a, S30, L0)
mods[[length(mods)+1]] <- psem(D4a, S40, L0)
mods[[length(mods)+1]] <- psem(D5a, S50, L0)
mods[[length(mods)+1]] <- psem(D6a, S60, L0)

## C: habitat specialization constrains leaf traits
mods[[length(mods)+1]] <- psem(S1a, D0a, L0)
mods[[length(mods)+1]] <- psem(S2a, D0a, L0)
mods[[length(mods)+1]] <- psem(S3a, D0a, L0)
mods[[length(mods)+1]] <- psem(S4a, D0a, L0)
mods[[length(mods)+1]] <- psem(S5a, D0a, L0)
mods[[length(mods)+1]] <- psem(S6a, D0a, L0)

## D: leaf traits determine productivity
mods[[length(mods)+1]] <- psem(L1, S10, D0a)
mods[[length(mods)+1]] <- psem(L2, S20, D0a)
mods[[length(mods)+1]] <- psem(L3, S30, D0a)
mods[[length(mods)+1]] <- psem(L4, S40, D0a)
mods[[length(mods)+1]] <- psem(LP, S50, D0a)
mods[[length(mods)+1]] <- psem(LQ, S60, D0a)

# Biological hypothesis code for single cause (null) models
BH <- c(BH, rep("iv",6), rep("viii",6), rep("ix",6), rep("v",6))
#BH <- c(BH, rep("viii",6), rep("ix",6))
length(BH); length(mods) #54

### plus plot type models (exact repeat of all the previous models)

### triangular models
## A: leaf traits determines both HS and CP, HS has an additional effect on CP
mods[[length(mods)+1]] <- psem(L8.1at, S1a)
mods[[length(mods)+1]] <- psem(L8.2at, D2a)
mods[[length(mods)+1]] <- psem(L8.3at, D3a)
mods[[length(mods)+1]] <- psem(L8.4at, D4a)
mods[[length(mods)+1]] <- psem(L8.5at, D5a)
mods[[length(mods)+1]] <- psem(L8.6at, D6a)

## B: HS determines leaf traits and CP, leaf traits have an additional effect on CP
mods[[length(mods)+1]] <- psem(L8.1at, S1a)
mods[[length(mods)+1]] <- psem(L8.2at, S2a)
mods[[length(mods)+1]] <- psem(L8.3at, S3a)
mods[[length(mods)+1]] <- psem(L8.4at, S4a)
mods[[length(mods)+1]] <- psem(L8.5at, S5a)
mods[[length(mods)+1]] <- psem(L8.6at, S6a)

# Biological hypothesis code for triangular models
BH <- c(BH, rep("i", 6), rep("ii", 6))

### chain models
## A: leaf traits determine specialization, which alters litter productivity
# i~55th
mods[[length(mods)+1]] <- psem(D1a, L5t)
mods[[length(mods)+1]] <- psem(D2a, L5t)
mods[[length(mods)+1]] <- psem(D3a, L5t)
mods[[length(mods)+1]] <- psem(D4a, L5t)
mods[[length(mods)+1]] <- psem(D5a, L5t)
mods[[length(mods)+1]] <- psem(D6a, L5t)

## B: specialization constrains leaf traits, which alter litter productivity
# i~61st
mods[[length(mods)+1]] <- psem(S1a, L1t)
mods[[length(mods)+1]] <- psem(S2a, L2t)
mods[[length(mods)+1]] <- psem(S3a, L3t)
mods[[length(mods)+1]] <- psem(S4a, L4t)
mods[[length(mods)+1]] <- psem(S5a, LPt)
mods[[length(mods)+1]] <- psem(S6a, LQt)

# Biological hypothesis code for chain models
BH <- c(BH, rep("vi", 6), rep("vii", 6))

### LAMBDA model
mods[[length(mods)+1]] <- psem(L8.1at, S10, D0a)
mods[[length(mods)+1]] <- psem(L8.2at, S20, D0a)
mods[[length(mods)+1]] <- psem(L8.3at, S30, D0a)
mods[[length(mods)+1]] <- psem(L8.4at, S40, D0a)
mods[[length(mods)+1]] <- psem(L8.5at, S50, D0a)
mods[[length(mods)+1]] <- psem(L8.6at, S60, D0a)

# Biological hypothesis code for lambda model
BH <- c(BH, rep("iii", 6))

### Single cause (null) models
## A: specialization determines litter production
mods[[length(mods)+1]] <- psem(L5t, D0a, S10)
mods[[length(mods)+1]] <- psem(L5t, D0a, S20)
mods[[length(mods)+1]] <- psem(L5t, D0a, S30)
mods[[length(mods)+1]] <- psem(L5t, D0a, S40)
mods[[length(mods)+1]] <- psem(L5t, D0a, S50)
mods[[length(mods)+1]] <- psem(L5t, D0a, S60)

## B: leaf traits determine distribution
mods[[length(mods)+1]] <- psem(D1a, S10, L0t)
mods[[length(mods)+1]] <- psem(D2a, S20, L0t)
mods[[length(mods)+1]] <- psem(D3a, S30, L0t)
mods[[length(mods)+1]] <- psem(D4a, S40, L0t)
mods[[length(mods)+1]] <- psem(D5a, S50, L0t)
mods[[length(mods)+1]] <- psem(D6a, S60, L0t)

## C: habitat specialization constrains leaf traits
mods[[length(mods)+1]] <- psem(S1a, D0a, L0t)
mods[[length(mods)+1]] <- psem(S2a, D0a, L0t)
mods[[length(mods)+1]] <- psem(S3a, D0a, L0t)
mods[[length(mods)+1]] <- psem(S4a, D0a, L0t)
mods[[length(mods)+1]] <- psem(S5a, D0a, L0t)
mods[[length(mods)+1]] <- psem(S6a, D0a, L0t)

## D: leaf traits determine productivity
mods[[length(mods)+1]] <- psem(L1t, S10, D0a)
mods[[length(mods)+1]] <- psem(L2t, S20, D0a)
mods[[length(mods)+1]] <- psem(L3t, S30, D0a)
mods[[length(mods)+1]] <- psem(L4t, S40, D0a)
mods[[length(mods)+1]] <- psem(LPt, S50, D0a)
mods[[length(mods)+1]] <- psem(LQt, S60, D0a)

# Biological hypothesis code for single cause (null) models
BH <- c(BH, rep("iv",6), rep("viii",6), rep("ix",6), rep("v",6))
#BH <- c(BH, rep("viii",6), rep("ix",6))

length(BH); length(mods) #108







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





###############################
#         APPENDIX B          #
# CANOPY PROD MODEL SELECTION #
###############################

library(MuMIn)
Lmods <- model.sel(L0,L0t,L1,L1t,L2,L2t,
  L3,L3t,L4,L4t,L5,L5t,L6,L6t,
  L7,L7t,L8.1a,L8.1at,L8.1b,L8.1bt,L8.1c,L8.1ct,
  L8.2a,L8.2at,L8.2b,L8.2bt,L8.2c,L8.2ct,L8.3a,L8.3at,
  L8.3b,L8.3bt,L8.3c,L8.3ct,L8.4a,L8.4at,L8.4b,L8.4bt,
  L8.4c,L8.4ct,L8.5a,L8.5at,L8.5b,L8.5bt,L8.5c,L8.5ct,
  L8.6a,L8.6at,L8.6b,L8.6bt,L8.6c,L8.6ct)

Lmods[Lmods$delta<4,]
summary(L2ti)
summary(L3ti)
summary(L8.2c)
summary(L8.3c)
summary(L8.2ct)
summary(L8.2a)
summary(L2)

# include models that dont contain basal area
summary(N0t <- lmer(loglitter ~ plotwet + (1|plot) + (logba|species), data=litter))
summary(N1t <- lmer(loglitter ~ plotwet + SLA + (1|plot) + (logba|species), data=litter))
summary(N2t <- lmer(loglitter ~ plotwet + LT + (1|plot) + (logba|species), data=litter))
summary(N3t <- lmer(loglitter ~ plotwet + LDMC + (1|plot) + (logba|species), data=litter))
summary(N4t <- lmer(loglitter ~ plotwet + CNratio + (1|plot) + (logba|species), data=litter))
summary(NPt <- lmer(loglitter ~ plotwet + PC1 + (1|plot) + (logba|species), data=litter))
summary(NQt <- lmer(loglitter ~ plotwet + PC2 + (1|plot) + (logba|species), data=litter))

summary(N0 <- lmer(loglitter ~ 1+ (1|plot) + (logba|species), data=litter))
summary(N1 <- lmer(loglitter ~ SLA + (1|plot) + (logba|species), data=litter))
summary(N2 <- lmer(loglitter ~ LT + (1|plot) + (logba|species), data=litter))
summary(N3 <- lmer(loglitter ~ LDMC + (1|plot) + (logba|species), data=litter))
summary(N4 <- lmer(loglitter ~ CNratio + (1|plot) + (logba|species), data=litter))
summary(NP <- lmer(loglitter ~ PC1 + (1|plot) + (logba|species), data=litter))
summary(NQ <- lmer(loglitter ~ PC2 + (1|plot) + (logba|species), data=litter))

summary(N7 <- lmer(loglitter ~ iv.wet + (1|plot) + (logba|species), data=litter))
summary(N8.1c <- lmer(loglitter ~ iv.wet + SLA + (1|plot) + (logba|species), data=litter))
summary(N8.2c <- lmer(loglitter ~ iv.wet + LT + (1|plot) + (logba|species), data=litter))
summary(N8.3c <- lmer(loglitter ~ iv.wet + LDMC + (1|plot) + (logba|species), data=litter))
summary(N8.4c <- lmer(loglitter ~ iv.wet + CNratio + (1|plot) + (logba|species), data=litter))
summary(N8.5c <- lmer(loglitter ~ iv.wet + PC1 + (1|plot) + (logba|species), data=litter))
summary(N8.6c <- lmer(loglitter ~ iv.wet + PC2 + (1|plot) + (logba|species), data=litter))

LIVmods <- model.sel(
  L0, L0t, L1, L1t, L2, L2t,
  L3, L3t, L4, L4t,
  L7, L8.1c, L8.2c, L8.3c,
  L8.4c, L8.5c, L8.6c,
  N0, N0t, N1, N1t, N2, N2t,
  N3, N3t, N4, N4t,
  N7, N8.1c, N8.2c, N8.3c,
  N8.4c, N8.5c, N8.6c
)
LIVmods

t2<-nrow(LIVmods)
top2<-vector("list",length=t2)

table <- data.frame(logba=rep(NA,t2),
  plotwet=NA,
  SLA=NA,
  LT=NA,
  LDMC=NA,
  CNratio=NA,
  iv.wet=NA,
  PC1=NA,
  PC2=NA)

for(i in 1:t2){
  coefmat <- summary(get.models(LIVmods, subset=i)[[1]])$coef
  for(j in names(table)){
    tryCatch({ 
    table[i,j] <- paste0(
      round(coefmat[j,"Estimate"], 2),
      " (?",
      round(coefmat[j,"Std. Error"], 2),
      ")")
    },error=function(e){})
} }

table <- 
cbind(table, LIVmods[,12:16])

for(i in 1:t2){
table$R2[i] <- 
 round(
  attr(
   r.squaredLR(
    get.models(LIVmods, subset=i)[[1]]),
  "adj.r.squared")
 ,3)
}
table[11:14] <- round(table[11:14],3)
table
#write.csv(table,"Table S3.csv")

#########################
# SUPPLEMENTARY FIGURES #
#########################

pal <- adjustcolor(rainbow(25), alpha.f=0.5)
litter$col <- as.numeric(as.factor(litter$species))

#jpeg("Fig. S2.jpg", res=300, height=8, width=15, units="in")
par(mar=c(5,5,2,2), mfrow=c(1,2))
plot(litter.production ~ ten.ba, data=litter, pch=16, col=pal[col], cex=4,
 ylab="Canopy production (log-transformed)", xlab="Plot summed basal area (log-transformed)",
 las=1, cex.lab=1.5, cex.axis=1.5)
sptab <- table(litter$species)
splist <- names(sptab)[which(sptab>2)]
for(i in 1:length(splist)){
 ltx <- tapply(litter$LT, litter$species, mean)[splist[i]]
 ivx <- tapply(litter$iv.wet, litter$species, mean)[splist[i]]
 bax <- with(subset(litter, species==splist[i]), seq(min(logba), max(logba), len=100))
 colx <- tapply(litter$col, litter$species, mean)[splist[i]]
 lines(
  exp(predict(L8.2c, newdata=data.frame(LT=ltx, logba=bax, species=splist[i], iv.wet=ivx), re.form= (~1|species)))
  ~ exp(bax, lwd=4, col=rainbow(25)[colx], lty=2)
}

plot(resid(L0) ~ iv.wet, data=litter, cex=4, cex.lab=1.5, cex.axis=1.5,
 xlab="Swamp indicator value", ylab="Residuals of basal area model", pch=16,
 col=pal[col], las=1)
abline(0, -1.8755, lty=2, lwd=5)
dev.off()

summary(L8.2c)
mean(litter$iv.wet)
mean(litter$LT)

brokenSIV <- cut(litter$iv.wet, 100)
levels(brokenSIV) <- seq(1,100,1)
col100 <- 101-as.numeric(as.character(brokenSIV))
colnamed <- tapply(col100, litter$species, mean)

library(stringr)
library(viridis)
#jpeg("Fig. S2 untransformed new.jpg", res=300, height=8, width=8, units="in")
par(mar=c(5,5,2,2))
plot(litter.production ~ ten.ba, data=litter, pch=16, cex=4, type="n", 
 ylab=expression(paste("Canopy production (g", " year"^"-1", " m"^"-2", ")")), 
 xlab=expression(paste("Plot summed basal area (",cm^2, ")")),
 las=1, cex.lab=1.5, cex.axis=1.5, col=adjustcolor(viridis(100)[col100], alpha.f=0.5))
sptab <- table(litter$species)
splist <- names(sptab)[which(sptab>0)]
spn <- names(sptab)[which(sptab>1)]
for(i in 1:length(splist)){
 ltx <- tapply(litter$LT, litter$species, mean)[splist[i]]
 ivx <- tapply(litter$iv.wet, litter$species, mean)[splist[i]]
 if(splist[i] %in% spn){
  bax <- with(subset(litter, species==splist[i]), 
   seq(min(logba), max(logba), len=100))
 } else {
  bax <- seq(0, max(litter$logba[which(litter$species==splist[i])]), len=100)
 }
 colx <- adjustcolor(viridis(100)[colnamed[match(splist[i], names(colnamed))]], alpha.f=0.7)
 lines(
  exp(predict(L8.2c, newdata=data.frame(LT=mean(litter$LT), logba=bax, species=splist[i], iv.wet=ivx), re.form= (~1|species)))
  ~ exp(bax), lwd=6, col=colx)
}
legend_image <- as.raster(matrix(viridis(20), ncol=1))
rasterImage(legend_image, 5800, 0, 6100, 50)
rect(5800, 0, 6100, 50)
text(rep(6150,6), seq(0,50,10), labels=format(seq(0,0.75,0.15), digits=2), adj=0, cex=1.2)
text(6150, 58, labels="SIV", adj=0.5, cex=1.5)
dev.off()


#jpeg("Fig. S2 untransformed new.jpg", res=300, height=8, width=8, units="in")
par(mar=c(5,5,2,2))
plot(litter.production ~ ten.ba, data=litter, pch=16, cex=2, xlim=c(0,8000),
 ylab=expression(paste("Canopy production (g", " year"^"-1", " m"^"-2", ")")), 
 xlab=expression(paste("Plot summed basal area (",cm^2, ")")),
 las=1, cex.lab=1.5, cex.axis=1.5, col=adjustcolor(viridis(100)[col100], alpha.f=0.4))
bax <- seq(0,6000,10)
ivx <- seq(0.01, 0.76, 0.25)
for(i in 1:4){
 lines(
  exp(predict(L7, newdata=data.frame(logba=bax, iv.wet=ivx[i]), re.form=NA))
  ~ exp(bax), lwd=6, col=viridis(4)[5-i])
}
legend_image <- as.raster(matrix(viridis(20), ncol=1))
rasterImage(legend_image, 6800, -5, 7100, 45)
rect(6800, -5, 7100, 45)
text(rep(7150,6), seq(-5,45,10), labels=format(seq(0,0.75,0.15), digits=2), adj=0, cex=1.2)
text(7150, 63, labels="Swamp\nIndicator Value", adj=0.5, cex=1.5)
dev.off()


















plot(resid(L0) ~ iv.wet, data=litter, cex=4, cex.lab=1.5, cex.axis=1.5,
 xlab="Swamp indicator value", ylab="Residuals of basal area model", pch=16,
 col=viridis(100)[col100], las=1)
abline(0, -1.8755, lty=2, lwd=5)

library(viridis)



litter[litter$iv.wet<0.1,]






