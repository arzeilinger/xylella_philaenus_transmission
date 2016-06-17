#### Daniele's Philaenus transmission study
#### 2015-11-19

# Preliminaries
rm(list = ls())
setwd("C:/Users/Adam/Documents/UC Berkeley post doc/Almeida lab/Daniele Philaenus study 2015")
my.packages <- c("lme4", "lmerTest", "lattice", "MASS", "tidyr", "dplyr", "nlme", "bbmle", "HH")
lapply(my.packages, require, character = TRUE)

# Load data
psdata <- read.csv("philaenus_data_qPCR.csv", header = TRUE)
str(psdata)


########################################################################
#### Transmission probability
#### Generalized linear mixed model, with source plant as random effect
## Model selection
# Model with only AAP
psModaap <- glmer(infected ~ log1p(xf.pop) + aap + (1|source.plant), data = psdata, family = "binomial")
# Model with only IAP
psModiap <- glmer(infected ~ log1p(xf.pop) + iap + (1|source.plant), data = psdata, family = "binomial")
# Model with both AAP and IAP
psMod <- glmer(infected ~ log1p(xf.pop) + aap + iap + (1|source.plant), data = psdata, family = "binomial")
# AIC model selection
AICctab(psMod, psModaap, psModiap, nobs = nrow(psdata), base = TRUE) # psMod is best
summary(psMod)


# AAP/IAP model with source.plant as fixed effect
# Set contrasts to "sum contrasts"
options(contrasts = c("contr.sum", "contr.poly"))
options("contrasts")
options(contrasts = c("contr.treatment", "contr.poly"))

psModFixed <- glm(infected ~ log1p(xf.pop) + aap + iap + source.plant, data = psdata, family = "binomial")
summary(psModFixed)


## Calculate proportion of non-infected test plants, for plotting
propInfectedAAP <- psdata[psdata$iap == 48,] %>% group_by(aap) %>% 
  summarise(n.infected = sum(infected, na.rm = TRUE), n = length(!is.na(infected)), prop.inf = n.infected/n)
propInfectedIAP <- psdata[psdata$aap == 48,] %>% group_by(iap) %>% 
  summarise(n.infected = sum(infected, na.rm = TRUE), n = length(!is.na(infected)), prop.inf = n.infected/n)

## Plotting results
# Remove NA (dead Ps)
psdata <- psdata[!is.na(psdata$infected),]
psdata <- psdata[!is.na(psdata$xf.pop),]
# Sample sizes across AAP and IAP treatments
table(psdata$aap, psdata$iap)
# Generate predicted values for both AAP and IAP
xv <- seq(0,50,length.out = nrow(psdata))
psModPlots <- glmer(infected ~ aap + iap + (1|source.plant), data = psdata, family = "binomial")
# Predictions for AAP
xaap <- data.frame(aap = xv,
                   iap = 48)
pred.aap <- predict(psModPlots, newdata = xaap, type = "response", re.form = NA)
# Predictions for IAP
xiap <- data.frame(aap = 48,
                   iap = xv)
pred.iap <- predict(psModPlots, newdata = xiap, type = "response", re.form = NA)


# Plot AAP and IAP results separately
pdf("Ps_duration_results_2plots_2016-04-06.pdf")#,
     #width = 76*2, height = 76, units = "mm", res = 600, compression = "lzw")
  par(mfrow=c(2,2))
  plot(propInfectedAAP$aap, propInfectedAAP$prop.inf, 
       col = "black", pch = 16,
       cex.axis = 1.1, cex.lab = 1.3, cex = 1.3,
       ylim = c(0,1), xlim = c(0,50),
       ylab = "Proportion of plants infected", xlab = "Acquisition access period (hr)")
  lines(smooth.spline(xv, pred.aap))
  plot(propInfectedIAP$iap, propInfectedIAP$prop.inf, 
       col = "black", pch = 16,
       cex.axis = 1.1, cex.lab = 1.3, cex = 1.3,
       ylim = c(0,1), xlim = c(0,50),
       ylab = "", 
       xlab = "Inoculation access period (hr)")
  lines(smooth.spline(xv, pred.iap))
dev.off()


# Exporting AAP and IAP figures as EPS


# # Plot AAP and IAP results together
# tiff("Ps_duration_results_plot_2016-01-14.tif")
#   plot(jitter(psdata$aap), psdata$infected, col = "darkgreen",
#        cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
#        ylim = c(0,1), xlim = c(0,50),
#        ylab = "Probability of transmission", xlab = "Duration (hrs)")
#   rug(jitter(psdata$iap[psdata$infected == 0]), col = "blue", lwd = 1)
#   rug(jitter(psdata$iap[psdata$infected == 1]), side = 3, col = "blue", lwd = 1)
#   lines(smooth.spline(xv, pred.aap), col = "darkgreen", lwd = 2, lty = 1) # P(transmission) for AAP
#   lines(smooth.spline(xv, pred.iap), col = "blue", lwd = 2, lty = 1) # P(transmission) for IAP
#   legend(x = 0, y = 0.9, legend = c("IAP", "AAP"), col = c("blue", "darkgreen"),
#          pch = c("|", "o"))
# dev.off()

#### Plotting Xf pop vs predicted transmission probability
psdata$predy <- predict(psMod, type = "response", re.form = NA)
psdata <- psdata[order(psdata$xf.pop),]
psdata$xf.popln <- log1p(psdata$xf.pop)
# Bin xf pop values and calculate the proportion of plants infected from each bin
xbin <- seq(0,8,by=1) # bin values for natural log xf pop
xfplot <- data.frame(xbin = xbin,
                     n = NA,
                     prop.inf = NA,
                     prop.inf.se = NA)
for(i in 1:length(xbin)){
  x.i <- xbin[i]
  x.ip1 <- xbin[i+1]
  infected.i <- psdata[psdata$xf.popln >= x.i & psdata$xf.popln < x.ip1,"infected"]
  n.i <- length(infected.i)
  prop.i <- sum(infected.i)/n.i
  prop.se.i <- sqrt(prop.i*(1-prop.i)/n.i)
  xfplot[i,"n"] <- n.i
  xfplot[i,"prop.inf"] <- prop.i
  xfplot[i,"prop.inf.se"] <- prop.se.i
}
xfplot <- xfplot[!is.nan(xfplot$prop.inf),]
# Make a new variable for the symbol expansion according to sample size of each bin
xfplot$nlog <- round(log10(xfplot$n))+1
xfplot$nlog[xfplot$nlog == 0] <- 0.5

pdf("Ps_xfpop_binned_plot_2016-04-06.pdf")
  plot(x = xfplot$xbin+0.5, y = xfplot$prop.inf, 
       cex.axis = 1.3, cex.lab = 1.3, cex = xfplot$nlog,
       ylim = c(0,1), xlim = c(0,8), pch = 16,
       ylab = "Proportion of plants infected", xlab = "Xylella population in vector (ln transformed)")
  lines(smooth.spline(log(psdata$xf.pop+1), psdata$predy, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()

tiff("Ps_xfpop_transmission_plot_2016-01-14.tif")
  plot(log(psdata$xf.pop+1), psdata$predy,
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       ylim = c(0,1), xlim = c(0,8),
       ylab = "Probability of transmission", xlab = "Xylella population in vector (log transformed)")
  lines(smooth.spline(log(psdata$xf.pop+1), psdata$predy, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()


##########################################################
#### Vector infection vs transmission
# Make a binary variable for vector infection
psdata$vector.infected <- psdata$xf.pop
# Turn vector.infected into a binary variable (1 = positive, 0 = negative)
psdata$vector.infected[psdata$vector.infected > 0] <- 1
table(psdata$vector.infected, psdata$infected)
chisq.test(psdata$vector.infected, psdata$infected)


#############################################################
#### Relating vector infection to AAP
## Function to run multiple models and compare them with AICc
modelComparison <- function(data = psdata){
  lmMod <- lme(log1p(xf.pop) ~ aap, random = ~ 1|source.plant, data = data)
  mmMod <- nlme(log1p(xf.pop) ~ a*aap/(b+aap), fixed = (a+b ~ 1), 
              random = a+b ~ 1|source.plant, data = data,
              start = c(a = 8, b = 4))
#   monoMod <- nlme(log1p(xf.pop) ~ a*(1-exp(-b*aap)), fixed = (a+b ~ 1), 
#               random = a+b ~ 1|source.plant, data = data,
#               start = c(a = 8, b = 4))
  aicComp <- AICctab(lmMod, mmMod, nobs = nrow(data), base = TRUE)
  resultsList <- list(aicComp, summary(lmMod), summary(mmMod), lmMod, mmMod)
  names(resultsList) <- c("AIC Comp", "Linear Model", "M-M Model", "Linear_Mod_Object", "MM_Mod_Object")
  return(resultsList)
}

resultsFull <- modelComparison(data = psdata)
results <- modelComparison(data = psdata[psdata$xf.pop > 0,])
# MM model is best

## Generating model fits from MM model

xv <- sort(runif(100, min = 0, max = 48))
# model fit with 0's
ahatFull <- as.numeric(fixef(resultsFull[["M-M Model"]]))[1]
bhatFull <- as.numeric(fixef(resultsFull[["M-M Model"]]))[2]
predMMFull <- ahatFull*xv/(bhatFull+xv)
# model fit without 0's
ahat <- as.numeric(fixef(results[["M-M Model"]]))[1]
bhat <- as.numeric(fixef(results[["M-M Model"]]))[2]
predMM <- ahat*xv/(bhat+xv)

pdf("xf_population_results_plot_2016-04-06.pdf")
  plot(jitter(psdata$aap, amount = 0), jitter(log1p(psdata$xf.pop), amount = 0),
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       xlim = c(0,50), ylim = c(0, 8),
       ylab = "Xylella population in vector (ln transformed)",
       xlab = "Acquisition access period (hr)")
  lines(xv, predMMFull, lwd = 2, lty = 1)
  lines(xv, predMM, lwd = 2, lty = 2)
  #abline(a = a, b = b, lwd = 2)
dev.off()


#### Evaluating fit of MM models
psdata$mmFullPred <- predict(resultsFull[["MM_Mod_Object"]])
summary(lm(log1p(xf.pop) ~ mmFullPred, data = psdata))
plot(x = psdata$mmFullPred, y = log1p(psdata$xf.pop))

infpsData <- psdata[psdata$xf.pop > 0,]
infpsData$mmPred <- predict(results[["MM_Mod_Object"]])
plot(x = infpsData$mmPred, y = log1p(infpsData$xf.pop))
summary(lm(log1p(xf.pop) ~ mmPred, data = infpsData))