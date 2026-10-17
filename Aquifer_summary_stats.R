########################################################
#### Microbes of the Morroa aquifer - summary stats ####
########################################################

rm(list=ls())
gc()

# Load packages
library(Hmsc)
library(abind)
library(stringr)
library(reshape2)
library(plyr)
library(readxl)
library(RgoogleMaps)

# Read and format taxonomy
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
write.csv(wide, "taxon_list.csv")

# Read assemblage data
dat = read.table("16S.txt", header=T)
cmat = dat[,-1]
cmat = data.frame(t(as.matrix(cmat)))

# Read metadata
mdat = read.table("meta_clean.txt", header=T)
summary(mdat)
names(mdat)

# Exclude blanks
rownames(cmat)
red = cmat[-c(27:30), ] # Blanks/control samples

drop = which(colSums(red>0)<1)
red = red[,-drop]
wide = wide[-drop,]

mdat = mdat[which(mdat$code_dna %in% row.names(red)),]
all(row.names(red) == mdat$code_dna)

# Taxomomic summary ####
nrow(wide)
sort(table(wide$d))
3632/3796 # 95.7% Bacteria

ar = wide[which(wide$d=="Archaea"),]
ba = wide[which(wide$d=="Bacteria"),]

sort(table(ba$p))

ps = wide[which(wide$p=="Pseudomonadota"),]
sort(table(ps$c))

sort(table(ar$p))

# Species richness summary
mean(rowSums(red>0))
median(rowSums(red>0))
sd(rowSums(red>0))
summary(rowSums(red>0))

tapply(rowSums(red>0), mdat$type, mean)
476.5714/123.4211
tapply(rowSums(red>0), mdat$type, sd)
tapply(rowSums(red>0), mdat$type, median)

# Indicator taxa ####

# Fecal indicators

# Enterobacteriaceae
wide[which(wide$f=="Enterobacteriacea"),]

sel = which(wide$f=="Enterobacteriacea")
reddat = data.frame(red[,sel], mdat$locality, mdat$depth_m, mdat$type)
head(reddat)

colSums(red[,sel]>1)
sum(reddat$X275>0)/sum(reddat$X275>-Inf)
tapply(reddat$X275, reddat$mdat.type, mean)

# Pseudomonas
wide[which(wide$g=="Pseudomonas"),]
unique(wide[which(wide$g=="Pseudomonas"),]$s)

sel = which(wide$g=="Pseudomonas")
length(sel)

reddat = data.frame(red[,sel], mdat$locality, mdat$depth_m, mdat$type, mdat$HDI)
head(reddat)

tapply(rowSums(red[,sel]>0), reddat$mdat.type, mean)
tapply(rowSums(red[,sel]), reddat$mdat.type, mean)

tapply(reddat[,1], reddat$mdat.type, mean)

comp = t(apply(reddat[,1:length(sel)], 2, function(x) tapply(x, reddat$mdat.type, mean)))
comp
plot(sqrt(comp[,1]), sqrt(comp[,2]), xlab="Abundance deep wells", ylab="Abundance piezometers", las=1)
lines(0:1000, 0:1000)

# Nitrosomonas (Nitrobacter not found)
wide[which(wide$g=="Nitrosomonas"),]
sel = which(wide$g=="Nitrosomonas")
length(sel)

reddat = data.frame(red[,sel], mdat$locality, mdat$depth_m, mdat$type, mdat$HDI)
reddat

# Hadarchaeota
sort(unique(ar$p))

wide[which(wide$p=="Hadarchaeota"),]
sel = which(wide$p=="Hadarchaeota")
length(sel)

reddat = data.frame(red[,sel], mdat$locality, mdat$depth_m, mdat$type, mdat$HDI)
reddat
tail(reddat, 10)
comp = t(apply(reddat[,1:length(sel)], 2, function(x) tapply(x, reddat$mdat.type, mean)))
comp

# Thermoproteota
sort(unique(ar$p))
wide[which(wide$p=="Thermoproteota"),]
sel = which(wide$p=="Thermoproteota")

reddat = data.frame(red[,sel], mdat$locality, mdat$depth_m, mdat$type, mdat$HDI)
reddat
tail(reddat, 10)

comp = t(apply(reddat[,1:length(sel)], 2, function(x) tapply(x, reddat$mdat.type, mean)))
comp
plot(sqrt(comp[,1]), sqrt(comp[,2]), xlab="Ind", ylab="Obs", las=1)
lines(0:1000, 0:1000)

sum(comp[,1]>0 & comp[,2]==0) # Ind but not Obs
sum(comp[,1]==0 & comp[,2]>0) # Obs but not Ind
sum(comp[,1]>0 & comp[,2]>0) # Both

# Basic figure
par(mfrow=c(1,1))
b = barplot(as.matrix(table(wide$d))/sum(as.matrix(table(wide$d))), xlim=c(0, 10), las=1)
axis(1, at=b, "Domain")
b = barplot(as.matrix(table(wide$p))/sum(as.matrix(table(wide$p))), xlim=c(0, 10), yaxt="n", space=2, add=TRUE)
axis(1, at=b, "Phylum")
b = barplot(as.matrix(table(wide$c))/sum(as.matrix(table(wide$c))), xlim=c(0, 10), yaxt="n", space=4, add=TRUE)
axis(1, at=b, "Class")
b = barplot(as.matrix(table(wide$o))/sum(as.matrix(table(wide$o))), xlim=c(0, 10), yaxt="n", space=6, add=TRUE)
axis(1, at=b, "Order")

# Complex figure
par(mfrow=c(1,1))
b = barplot(as.matrix(table(wide$d))/sum(as.matrix(table(wide$d))), xlim=c(0, 10), las=1)
axis(1, at=b, "Domain")

b = barplot(as.matrix(table(ba$p))/sum(as.matrix(table(ba$p))), xlim=c(0, 10), yaxt="n", space=2, add=TRUE)
axis(1, at=b, "Bacterial Phylum")

b = barplot(as.matrix(table(ps$c))/sum(as.matrix(table(ps$c))), xlim=c(0, 10), yaxt="n", space=4, add=TRUE)
axis(1, at=b, "Pseudomonodota Class")

