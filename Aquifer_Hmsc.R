########################################################
#### Microbes of the Morroa aquifer - Hmsc analysis ####
########################################################

rm(list=ls())
gc()

# Load packages
library(readxl)
library(stringr)
library(plyr)
library(reshape2)
library(Hmsc)
library(abind)
library(fBasics)

# Read data
taxonomy = data.frame(read_excel("taxonomy.xlsx"))
out = strsplit(taxonomy$taxonomy, split=";")

for(i in 1:length(out)){
  out[[i]] = gsub(" ", "", out[[i]])
  out[[i]] = data.frame(ID=rep(i, length(out[[i]])), level = substr(out[[i]], 1,1), taxon = substr(out[[i]], 4, 20))
}

out = rbind.fill(out)
out$level = factor(out$level, levels=c("d", "p", "c", "o", "f", "g", "s"))
wide = dcast(out, ID~level)
head(wide)

# Community data
dat <- read.table("16S.txt", header=T)
cmat = dat[,-1]
cmat = data.frame(t(as.matrix(cmat)))

mdat = read.table("meta_clean.txt", header=T)
summary(mdat)
names(mdat)

# Exclude blanks etc. and singletons
rownames(cmat)
red = cmat[-c(27:30), ]

drop = which(colSums(red>0)<2)
red = red[,-drop]
wide = wide[-drop,]

# Hmsc setup full data ####

# Construct Y
Y = red
dim(Y)

par(mfrow=c(1,2), mar=c(4,4,2,2))
hist(colSums(Y>0), main="Sites per species", xlab="Sites")
hist(rowSums(Y>0), main="Species per site", xlab="Species")

# Construct XData
mdat = mdat[which(mdat$code_dna %in% row.names(Y)),]
all(row.names(Y) == mdat$code_dna)

XData = subset(mdat, select=c("depth_m", "type", "C_L_C", "HDI", "long", "lat"))
XData$type = as.factor(XData$type)
XData$C_L_C = as.factor(XData$C_L_C)
XData$effort = rowSums(Y)
summary(XData)

# Construct TrData
names(wide)
TrData = subset(wide, select=c("d", "p"))
rownames(TrData) = colnames(Y)
head(TrData)

# Construct dfPi
dfPi = data.frame(SU=as.factor(row.names(red)))
all(dfPi$SU == row.names(Y))

# Set random levels
rL1 = HmscRandomLevel(units = unique(dfPi$SU))
rL1

# Set model formulae
names(XData)
XFormula = ~ depth_m + HDI + type
TrFormula= ~ p

# Setup HMSC models
m1 = Hmsc(Y = as.matrix((Y>0)*1), 
          XData = XData,  XFormula = XFormula,
          TrData = TrData, TrFormula = TrFormula,
          distr = "probit", 
          studyDesign = dfPi, 
          ranLevels = list(SU=rL1))

# Sample MCMC
thin = 300
samples = 250
adaptNf = ceiling(0.4*samples*thin)
transient = ceiling(0.5*samples*thin)
nChains = 2

a = Sys.time()
m1 = sampleMcmc(m1, samples = samples, thin = thin,
                adaptNf = rep(adaptNf, 1),
                transient = transient,
                nChains = nChains, nParallel = 2, updater=list(GammaEta=FALSE))
Sys.time()-a

# 764 species, thin 300, 2 chains (parallel): 8.2 hours

save(m1, file="PAmod/fitted_mods/mPA_thin300_2chains.Rdata")

#### Presence-absence model postprocessing
load(file="PAmod/fitted_mods/mPA_thin300_2chains.Rdata")

# Assess sampling performance
mpost = convertToCodaObject(m1)

esBeta = effectiveSize(mpost$Beta)
summary(esBeta)

psrf = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
summary(psrf)

# Plot posterior trace
pdf("PAmod/betaPostAll.pdf")
plot(mpost$Beta[,1:200])
dev.off()

# Explanatory power
predY = computePredictedValues(m1, expected=T)
MF = evaluateModelFit(hM=m1, predY = predY)

mean(MF$TjurR2, na.rm=T)
range(MF$TjurR2, na.rm=T)

mean(MF$AUC, na.rm=T)
range(MF$AUC, na.rm=T)

par(mfrow=c(1,1))
plot(colSums(m1$Y)/nrow(m1$Y), MF$TjurR2, pch=16, las=1,
     xlab="Prevalence", ylab="Coefficient of Discrimination (Tjur's D)")

predYm = apply(predY, 1:2, mean)

plot(colSums(predYm,na.rm=T), colSums(m1$Y))
lines(-1000:1000, -1000:1000)
cor(colSums(predYm), colSums(m1$Y))^2

plot(rowSums(predYm), rowSums(m1$Y), las=1,
     xlab="Predicted species richness",
     ylab="Observed species richness")
lines(-1000:1000, -1000:1000)
cor(rowSums(predYm), rowSums(m1$Y))^2

# Variance partitioning ####
cbind(m1$covNames)
group = c(rep(1,2), 2, 3)
cbind(m1$covNames, group)

groupnames = c("Depth", "HDI", "WellType")
VP = computeVariancePartitioning(hM = m1, group = group, groupnames = groupnames)
rowMeans(VP$vals)*100
apply(VP$vals, 1, range)*100
VP$R2T #Phylum effects

VP$vals = VP$vals[,rev(order(colSums(VP$vals[1:3,])))]

pdf("PAmod/VarPartAll.pdf", height=3, width=30, family="Times")
par(mar=c(2,4,2,12))
plotVariancePartitioning(m1, VP)
dev.off()

# Beta parameters
pbeta = getPostEstimate(m1, "Beta")
m1$covNames

HDIp = wide[which(pbeta$support[3,]>.9),] # Positive response to HDI
dim(HDIp)

HDIn = wide[which(pbeta$support[3,]<.1),] # Positive response to HDI
dim(HDIn)

x11()
plotBeta(m1, post=pbeta, spNamesNumbers = c(F, F))

# Gamma parameters
pgamma = getPostEstimate(m1, "Gamma")

x11()
plotGamma(m1, post=pgamma)

# Hmsc setup piezometers ####

# Construct Y
Y = red
dim(Y)

par(mfrow=c(1,2), mar=c(4,4,2,2))
hist(colSums(Y>0), main="Sites per species", xlab="Sites")
hist(rowSums(Y>0), main="Species per site", xlab="Species")

# Construct XData
mdat = mdat[which(mdat$code_dna %in% row.names(Y)),]
all(row.names(Y) == mdat$code_dna)

XData = subset(mdat, select=c("depth_m", "type", "C_L_C", "HDI", "long", "lat"))
XData$type = as.factor(XData$type)
XData$C_L_C = as.factor(XData$C_L_C)
XData$effort = rowSums(Y)
summary(XData)

# Construct TrData
names(wide)
TrData = subset(wide, select=c("d", "p"))
rownames(TrData) = colnames(Y)
head(TrData)

# Construct dfPi
dfPi = data.frame(SU=as.factor(row.names(red)))
all(dfPi$SU == row.names(Y))

# Set random levels
rL1 = HmscRandomLevel(units = unique(dfPi$SU))
rL1

# Set model formulae
names(XData)
XFormula = ~ depth_m + HDI

TrFormula= ~ p

# Drop deep wells
drop = which(XData$type=="Ind")
Y = Y[-drop,]
drop2 = which(colSums(Y)==0)
Y = Y[,-drop2]

TrData = TrData[-drop2,]
XData = XData[-drop,]

dfPi = data.frame(SU=as.factor(row.names(red[-drop,])))
all(dfPi$SU == row.names(Y))
rL1 = HmscRandomLevel(units = unique(dfPi$SU))
rL1

# Setup HMSC models
m1 = Hmsc(Y = as.matrix((Y>0)*1), 
          XData = XData,  XFormula = XFormula,
          TrData = TrData, TrFormula = TrFormula,
          distr = "probit", 
          studyDesign = dfPi, 
          ranLevels = list(SU=rL1))
m1

# Sample MCMC
thin = 500
samples = 250
adaptNf = ceiling(0.4*samples*thin)
transient = ceiling(0.5*samples*thin)
nChains = 2

a = Sys.time()
m1 = sampleMcmc(m1, samples = samples, thin = thin,
                adaptNf = rep(adaptNf, 1),
                transient = transient,
                nChains = nChains, nParallel = 2, updater=list(GammaEta=FALSE))
Sys.time()-a

save(m1, file="PAmod/fitted_mods/mPA_piezo_thin500_2chains.Rdata")

# Presence-absence model postprocessing Piezometers ####
load(file="PAmod/fitted_mods/mPA_piezo_thin500_2chains.Rdata")

# Assess sampling performance
mpost = convertToCodaObject(m1)

esBeta = effectiveSize(mpost$Beta)
summary(esBeta)

psrf = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
summary(psrf)

# Plot posterior trace
pdf("PAmod/betaPostPiezo.pdf")
plot(mpost$Beta[,1:200])
dev.off()

# Explanatory power
predY = computePredictedValues(m1, expected=T)
MF = evaluateModelFit(hM=m1, predY = predY)

mean(MF$TjurR2, na.rm=T)
range(MF$TjurR2, na.rm=T)

mean(MF$AUC, na.rm=T)
range(MF$AUC, na.rm=T)

pdf("TjurPlot.pdf", height=4.1, width=4.5, fam="Times")
par(mfrow=c(1,1))
plot(colSums(m1$Y)/nrow(m1$Y), MF$TjurR2, pch=16, las=1,
     xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
mtext("Prevalence", 1, line=2.5)
mtext("Coefficient of Discrimination (D)", 2, line=2.5)
dev.off()

predYm = apply(predY, 1:2, mean)

plot(colSums(predYm,na.rm=T), colSums(m1$Y))
lines(-1000:1000, -1000:1000)
cor(colSums(predYm), colSums(m1$Y))^2

plot(rowSums(predYm), rowSums(m1$Y), las=1,
     xlab="Predicted species richness",
     ylab="Observed species richness")
lines(-1000:1000, -1000:1000)
cor(rowSums(predYm), rowSums(m1$Y))^2

#### Variance partitioning ####
cbind(m1$covNames)
group = c(rep(1,2), 2)
cbind(m1$covNames, group)

groupnames = c("Depth", "HDI")
VP = computeVariancePartitioning(hM = m1, group = group, groupnames = groupnames)
str(VP)
rowMeans(VP$vals)*100
apply(VP$vals, 1, range)*100
VP$R2T #Phylum effects

VP$vals = VP$vals[,rev(order(colSums(VP$vals[1:3,])))]

pdf("PAmod/VarPartPiezo.pdf", height=3, width=30, family="Times")
par(mar=c(2,4,2,12))
plotVariancePartitioning(m1, VP)
dev.off()

# Beta parameters
pbeta = getPostEstimate(m1, "Beta")
str(pbeta)

m1$covNames

HDIp = wide[which(pbeta$support[3,]>.9),] # Positive response to HDI
dim(HDIp)

sort(table(HDIp$d))
ba = HDIp[which(HDIp$d=="Bacteria"),]
sort(table(ba$p))

HDIn = wide[which(pbeta$support[3,]<.1),] # Negative response to HDI
dim(HDIn)

sort(table(HDIn$d))
ba = HDIn[which(HDIn$d=="Bacteria"),]

sort(table(ba$p))

ps = ba[which(ba$p=="Pseudomonadota"),]

sort(table(ps$f))

write.csv(HDIn, file="Piezometers_Negative_HDI.csv")
write.csv(HDIp, file="Piezometers_Positive_HDI.csv")

x11()
plotBeta(m1, post=pbeta, spNamesNumbers = c(F, F))

# Gamma parameters
pgamma = getPostEstimate(m1, "Gamma")

x11()
plotGamma(m1, post=pgamma)

# Gradient plots ####
m1$covNames

Gradient = constructGradient(m1, focalVariable = "HDI", 
                             non.focalVariables = list(
                              depth_m = list(1)), 
                             ngrid=300)
head(Gradient$XDataNew, 10)

predY = predict(m1, Gradient=Gradient, expected=TRUE, predictEtaMean=FALSE)

plotGradient(m1, Gradient, predY, measure="S", showData = T)
plotGradient(m1, Gradient, predY, measure="Y", index=30, showData = T)

pdf("PAmod/HDIGradients.pdf", height=3.4, width=6.6, family="Times", pointsize=11)
par(mfrow=c(1,2), mar=c(4,4,2,1), oma=c(0,0,0,0))

xx = Gradient$XDataNew[,1]

# Per taxon
cols = 1:1000

pd = data.frame(Intercept=1, Gradient$XDataNew[,2], Gradient$XDataNew[,1])
betas = pbeta$mean[,1]
preds = pnorm(as.matrix(pd) %*% as.matrix(betas))

plot(Gradient$XDataNew$HDI, preds, ylim=c(0, 0.4), col="white", las=1,
     ylab="", xlab="")
mtext("Human Disturbance", 1, line=2.5)
mtext("Occurrence probability", 2, line=2.5, xpd=T)

sel = which(pbeta$support[3,]>.9 | pbeta$support[3,]<0.1)

for(i in sel){
  betas = pbeta$mean[,i]
  preds = pnorm(as.matrix(pd) %*% as.matrix(betas))
  lines(Gradient$XDataNew$HDI, preds, col=cols[i], lwd=2)
}

cols = topo.colors(4, alpha=.5)

predS = lapply(predY, function(x) rowSums(x))
pred = apply(abind(predS, along = 2), c(1), quantile, 
             prob = c(.025, .5, .975), na.rm = TRUE)
range(pred[2,])
plot(m1$XData$HDI, rowSums(m1$Y), pch=16, lwd = 2, col="grey", las=1,
     ylab="", xlab="")
mtext("Human Disturbance", 1, line=2.5)
mtext("Taxon richness", 2, line=2.5, xpd=T)

lines(xx, pred[2, ], lwd=2)
polygon(c(xx, rev(xx)), c(pred[1, ], rev(pred[3, ])), 
        col = cols[3], border = FALSE)

dev.off()

# Phylum plot
pdf("PAmod/PhylumGradients.pdf", height=3.4, width=5, family="Times", pointsize=11)

phyla = levels(as.factor(m1$TrData$p))[c(which(table(m1$TrData$p)>2))]
phyla
sort(table(m1$TrData$p))

predS = lapply(predY, function(x) rowSums(x[,which(m1$TrData$p==phyla[3])]))
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
range(egpred[2,]) # Bacillota

predS = lapply(predY, function(x) rowSums(x[,which(m1$TrData$p==phyla[8])]))
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
range(egpred[2,]) # Pseudomonadota

predS = lapply(predY, function(x) rowSums(x[,which(m1$TrData$p==phyla[3])]))
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
range(egpred[2,]) # Bacillota

#cols = topo.colors(10)
cols = divPalette(10, "Spectral")

predS = lapply(predY, function(x) rowSums(x[,which(m1$TrData$p==phyla[1])]))
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
par(mar=c(4,4,2,11))
plot(xx, log10(egpred[2, ]), type="l",lwd = 2, las=1, ylim=c(-2,2),
     ylab="", xlab="", col=cols[i], yaxt="n")
mtext("Human disturbance", 1, line=2.5)
mtext("Taxon richness", 2, line=2.5)
axis(2, seq(-2, 2, 1), round(10^seq(-2, 2, 1), 2), las=1)
#polygon(c(xx, rev(xx)), c(egpred[1, ], rev(egpred[3, ])), 
#        col = cols[1], border = FALSE)

for(i in 2:length(phyla)){
  predS = lapply(predY, function(x) rowSums(x[,which(m1$TrData$p==phyla[i])]))
  elpred = apply(abind(predS, along = 2), c(1), quantile, 
                 prob = c(.025, .5, .975), na.rm = TRUE)
  lines(xx, log10(elpred[2, ]), lwd = 2, col=cols[i])
  #polygon(c(xx, rev(xx)), c(elpred[1, ], rev(elpred[3, ])), 
   #       col = cols[i], border = FALSE)
}

par(xpd=T)
legend(5.2, 2, phyla, cols, bty="n")

dev.off()
