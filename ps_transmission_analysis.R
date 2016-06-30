#### Philaenus transmission study
#### 2016-06-17
#### Code written by Adam R. Zeilinger

#### Preliminaries
# Clear workspace
rm(list = ls())

# Load required packages
my.packages <- c("lme4", "lmerTest", "lattice", "MASS", "tidyr", 
                 "dplyr", "nlme", "bbmle", "HH", 'RCurl', 'foreign')
lapply(my.packages, require, character = TRUE)

# Load data from Github repository
url <- "https://raw.githubusercontent.com/arzeilinger/xylella_philaenus_transmission/master/SupplMat.csv"
psdata <- getURL(url) %>% textConnection() %>% read.csv(., header = TRUE)
str(psdata)


########################################################################
#### Transmission probability vs. AAP, IAP, and Xylella populations in vectors
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
pdf("Figure_1_Ps_duration_results_2plots.pdf")#,
     #width = 76*2, height = 76, units = "mm", res = 600, compression = "lzw")
#trellis.device()
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
#export.eps("Figure_1_Ps_duration_results_2plots.eps", width = 10, height = 7)
dev.off()

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
# Shift values to median of binned intervals
xfplot$nlog[xfplot$nlog == 0] <- 0.5

#tiff("Figure_2_Ps_xfpop_binned_plot.tif")
trellis.device()
  plot(x = xfplot$xbin+0.5, y = xfplot$prop.inf, 
       cex.axis = 1.3, cex.lab = 1.3, cex = xfplot$nlog,
       ylim = c(0,1), xlim = c(0,8), pch = 16,
       ylab = "Proportion of plants infected", xlab = "Xylella population in vector (ln transformed)")
  lines(smooth.spline(log(psdata$xf.pop+1), psdata$predy, nknots = 4, tol = 1e-6), lwd = 2)
export.eps("Figure_2_Ps_xfpop_binned_plot.eps")
#dev.off()


##############################################################################
#### Contingency table and Chi-square test of Vector and test plant infection status
# Make a binary variable for vector infection
psdata$vector.infected <- psdata$xf.pop
# Turn vector.infected into a binary variable (1 = positive, 0 = negative)
psdata$vector.infected[psdata$vector.infected > 0] <- 1

# Contingency table of vector and plant infection status
# 1 = infected/PCR-positive; 0 = non-infected/PCR-negative
table(psdata$vector.infected, psdata$infected, dnn = c("vector", "plant"))
# Just vector infection status
table(psdata$vector.infected)
# Just test plant infection status
table(psdata$infected)

# Chi-square test of contingency table
chisq.test(psdata$vector.infected, psdata$infected)


#############################################################
#### Relating vector infection to AAP
## Function to run multiple models and compare them with AICc
modelComparison <- function(data = psdata){
  lmMod <- lme(log1p(xf.pop) ~ aap, random = ~ 1|source.plant, data = data)
  mmMod <- nlme(log1p(xf.pop) ~ a*aap/(b+aap), fixed = (a+b ~ 1), 
              random = a+b ~ 1|source.plant, data = data,
              start = c(a = 8, b = 4))
  aicComp <- AICctab(lmMod, mmMod, nobs = nrow(data), base = TRUE)
  resultsList <- list(aicComp, summary(lmMod), summary(mmMod), lmMod, mmMod)
  names(resultsList) <- c("AIC Comp", "Linear Model", "M-M Model", "Linear_Mod_Object", "MM_Mod_Object")
  return(resultsList)
}

# Run multi-model comparison for entire data set (PCR-positive and PCR-negative vectors)
resultsFull <- modelComparison(data = psdata)

# Run multi-model comparison for data set of only PCR-positive vectors
results <- modelComparison(data = psdata[psdata$xf.pop > 0,])
# MM model is best in bot cases

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

#tiff("Figure_3_xf_population_results_plot.tif")
trellis.device()
  plot(jitter(psdata$aap, amount = 0), jitter(log1p(psdata$xf.pop), amount = 0),
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       xlim = c(0,50), ylim = c(0, 8),
       ylab = "Xylella population in vector (ln transformed)",
       xlab = "Time after beginning of aquisition access period (hr)")
  lines(xv, predMMFull, lwd = 2, lty = 1)
  lines(xv, predMM, lwd = 2, lty = 2)
export.eps("Figure_3_xf_population_results_plot.eps")
#dev.off()


#### Evaluating fit of MM models
# Full data set
psdata$mmFullPred <- predict(resultsFull[["MM_Mod_Object"]]) # Generate Xylella population predictions
summary(lm(log1p(xf.pop) ~ mmFullPred, data = psdata))
plot(x = psdata$mmFullPred, y = log1p(psdata$xf.pop),
     xlab = "Predicted Xylella populations",
     ylab = "Observed Xylella populations")

# Only PCR-positive ("infected") vectors
infpsData <- psdata[psdata$xf.pop > 0,]
infpsData$mmPred <- predict(results[["MM_Mod_Object"]]) # Generate Xylella population predictions
summary(lm(log1p(xf.pop) ~ mmPred, data = infpsData))
plot(x = infpsData$mmPred, y = log1p(infpsData$xf.pop),
     xlab = "Predicted Xylella populations",
     ylab = "Observed Xylella populations")


#################################################################################
#### Total latent time vs. xf populations
#### total latent time = AAP + IAP
psdata$latent.time <- psdata$aap + psdata$iap

modelComparisonLatent <- function(data = psdata){
  lmMod <- lme(log1p(xf.pop) ~ latent.time, random = ~ 1|source.plant, data = data)
  mmMod <- nlme(log1p(xf.pop) ~ a*latent.time/(b+latent.time), fixed = (a+b ~ 1), 
              random = a+b ~ 1|source.plant, data = data,
              start = c(a = 8, b = 4))
  aicComp <- AICctab(lmMod, mmMod, nobs = nrow(data), base = TRUE)
  resultsList <- list(aicComp, summary(lmMod), summary(mmMod), lmMod, mmMod)
  names(resultsList) <- c("AIC Comp", "Linear Model", "M-M Model", "Linear_Mod_Object", "MM_Mod_Object")
  return(resultsList)
}

# Run multi-model comparison for entire data set (PCR-positive and PCR-negative vectors)
resultsFull <- modelComparisonLatent(data = psdata)

# Run multi-model comparison for data set of only PCR-positive vectors
results <- modelComparisonLatent(data = psdata[psdata$xf.pop > 0,])
# MM model is best in bot cases

## Generating model fits from MM model

xv <- sort(runif(100, min = 0, max = max(psdata$latent.time, na.rm = TRUE)))
# model fit with 0's
ahatFull <- as.numeric(fixef(resultsFull[["M-M Model"]]))[1]
bhatFull <- as.numeric(fixef(resultsFull[["M-M Model"]]))[2]
predMMFull <- ahatFull*xv/(bhatFull+xv)
# model fit without 0's
ahat <- as.numeric(fixef(results[["M-M Model"]]))[1]
bhat <- as.numeric(fixef(results[["M-M Model"]]))[2]
predMM <- ahat*xv/(bhat+xv)

#tiff("xf_population_vs_latent_period_plot.tif")
trellis.device()
  plot(jitter(psdata$latent.time, amount = 0), jitter(log1p(psdata$xf.pop), amount = 0),
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       xlim = c(40, 100), ylim = c(0, 8),
       ylab = "Xylella population in vector (ln transformed)",
       xlab = "Time after beginning of aquisition access period (hr)")
  lines(xv, predMMFull, lwd = 2, lty = 1)
  lines(xv, predMM, lwd = 2, lty = 2)
export.eps("Figure_3_xf_population_latent_time_plot.eps")
#dev.off()


#### Evaluating fit of MM models
# Full data set
psdata$mmFullPred <- predict(resultsFull[["MM_Mod_Object"]]) # Generate Xylella population predictions
summary(lm(log1p(xf.pop) ~ mmFullPred, data = psdata))
plot(x = psdata$mmFullPred, y = log1p(psdata$xf.pop),
     xlab = "Predicted Xylella populations",
     ylab = "Observed Xylella populations")

# Only PCR-positive ("infected") vectors
infpsData <- psdata[psdata$xf.pop > 0,]
infpsData$mmPred <- predict(results[["MM_Mod_Object"]]) # Generate Xylella population predictions
summary(lm(log1p(xf.pop) ~ mmPred, data = infpsData))
plot(x = infpsData$mmPred, y = log1p(infpsData$xf.pop),
     xlab = "Predicted Xylella populations",
     ylab = "Observed Xylella populations")
