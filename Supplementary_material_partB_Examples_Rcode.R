
# The three example can be conveniently reproduced in R, relying only on
# the add-on packages mvtnorm, multcomp, mratios, and ggplot2 as well as 
# packages MASS to load the data of case study 1 and drc to load the data of
# case study 3.

library(multcomp)
library(mratios)
library(ggplot2)


###############################################
# Example 4.1: All pairwise comparisons in    #
# a linear model with baseline as covariate   #
###############################################

# Loading the data set
library(MASS)
data(anorexia)
AN <- anorexia
fitAN <- lm(Postwt ~ Treat-1 + Treat:Prewt, data=AN)

# Set of model predictions at covariate values 
VX <- matrix(seq(70,95,5), ncol=1)
ITx <- matrix(rep(1,length(VX)),ncol=1)
XM <- cbind(ITx, VX)

# Matrices for all pairwise comparisons
AMap <- contrMat(n=table(AN$Treat), type="Tukey")
CMap <- kronecker(XM, AMap)

# Compute the intervals 
COMPap <- glht(fitAN, linfct=CMap)
SCIap <- confint(COMPap)

# Plot corresponding to Figure 2)
dtrtap  <-  factor(rep(rownames(AMap), times=length(VX)), levels=rownames(AMap))
dvxap  <-  rep(VX, each=nrow(AMap))
dap <- data.frame(SCIap$confint, comp=dtrtap, prewt=dvxap)

ggplot(dap, aes(x=prewt, y=Estimate, ymin=lwr, ymax=upr)) + 
 geom_errorbar(width=1) + geom_line() + geom_point(shape=15) +
 facet_grid(.~comp) + xlab("Preweight") + 
 ylab("Difference in expected postweight") + geom_hline(yintercept=0)


###############################################
# Example 4.2: Comparisons to control with    #
# interaction to a quadratic regression term  #
###############################################

# The initial lines of code create a simulated data set that resembles the
# structure of the data set presented in Milliken and Johnson (2002). However,
# the data values are simulated, thus the following code will NOT exactly
# reproduce the estimates and confidence limits shown in the manuscript for
# this example. Reference: Milliken GA and Johnson DE (2002): Analysis of 
# Messy Data, Vol.II: Analysis of Covariance. Chapman & Hall/HRC.

rm(list=ls())

PC  <-  structure(list(yield = c(1.9, 4.5, 3, 2.2, 3.7, 4.3, 3.1, 5.6, 6.1, 
 5.5, 6, 5.8, 7.3, 7.7, 6.4, 7.8, 7.2, 6.3, 6.2, 5.9, 4.6, 5.1, 4.5, 2.8, 8.6,
 10, 8, 9, 7.1, 9.2, 9.4, 7.4, 7.6, 6.7, 9.1, 8.8), x = c(1.2, 3.3, 1.7, 2,
 2.2, 3.9, 1.9,  8.2, 6.9, 8.1, 7.7, 6.6, 4.6, 4.6, 2.7, 4.8, 4.2, 2.6, 3, 7.7,
 8.7, 8.8, 8.8, 9.8, 2.7, 4.3, 1.5, 2.5, 0.8, 3.5, 2.7, 8, 8.3, 8.6, 6.4, 6.5),
 additive = factor(rep(c("Control", "S1", "S2"), c(12,12,12)))),
 .Names = c("yield", "x", "additive"), class = "data.frame", row.names = c(NA, -36L))

# Model selected by Milliken and Johnson:
fitPC <- lm(yield ~ additive-1 + x + additive:I(x^2), data=PC)
coef(fitPC)


# Differences to control

# Matrices for differences to control
CTm1 <- contrMat(n=table(PC$additive), type="Dunnett")

# Grid of covariate values
VX <- matrix(seq(0,10,1), ncol=1)
Im1 <- matrix(rep(1, nrow(CTm1)), ncol=1)
ITx <- matrix(rep(1,length(VX)),ncol=1)
C1 <- kronecker(ITx, CTm1)
C2 <- kronecker(VX, Im1)*0
C3 <- kronecker(VX^2, CTm1)
CMm1 <- cbind(C1,C2,C3)
CMm1

# Compute the confidence intervals
COMPm1 <- glht(fitPC, linfct=CMm1)
SCIm1 <- confint(COMPm1)

# convert to a data.frame for plotting
dtrtm1  <-  factor(rep(rownames(CTm1), times=length(VX)), levels=rownames(CTm1))
dvxm1  <-  rep(VX, each=nrow(CTm1))
dm1 <- data.frame(SCIm1$confint, comp=dtrtm1, x=dvxm1)

# Plot corresponding to Figure 5)
ggplot(dm1, aes(x=x, y=Estimate, ymin=lwr, ymax=upr)) +
 geom_errorbar(width=0.3) + geom_line() + geom_point(shape=15) + 
 facet_grid(.~comp) + xlab("Substance x") + ylab("Difference in expected yield") +
 geom_hline(yintercept=0)


# Ratios to control

# Matrices for confidence intervals for ratios to control (relying on mratios)
ABm1 <- contrMatRatio(n=table(PC$additive), type="Dunnett")
Am1 <- ABm1$numC
Bm1 <- ABm1$denC
VX <- matrix(seq(1,9,1), ncol=1)
Im1 <- matrix(rep(1, nrow(Am1)), ncol=1)
ITx <- matrix(rep(1,length(VX)),ncol=1)
C1 <- kronecker(ITx, Am1)
C2 <- kronecker(VX, Im1)
C3 <- kronecker(VX^2, Am1)
CMm1 <- cbind(C1,C2,C3)
D1 <- kronecker(ITx, Bm1)
D2 <- kronecker(VX, Im1)
D3 <- kronecker(VX^2, Bm1)
DMm1 <- cbind(D1,D2,D3)
modelest  <-  coef(fitPC)
modelvc  <-  vcov(fitPC)

# Compute confidence intervals
RSCIm1 <- gsci.ratio(est=modelest, vcmat=modelvc,
 Num.Contrast=CMm1, Den.Contrast=DMm1, degfree = 29)

RSCIm1

# Plot corresponding to Figure 6)
dtrtm1  <-  factor(rep(rownames(Am1), times=length(VX)), levels=rownames(Am1))
dvxm1  <-  rep(VX, each=nrow(Am1))
dm1 <- data.frame(RSCIm1$estimate, RSCIm1$conf.int, comp=dtrtm1, x=dvxm1)

ggplot(dm1, aes(x=x, y=estimate, ymin=lower, ymax=upper, width=0.5)) +
 geom_errorbar(width=0.3) + geom_line() + geom_point(shape=15) +
 facet_grid(.~comp) + xlab("Substance x") + ylab("Ratio of expected yield") +
 geom_hline(yintercept=1) + scale_y_log10(breaks=c(0.5, 0.8, 1, 1.2, 1.5, 2, 5, 10))


######################################################
# Example 4.3: All pairwise comparisons in a         #
# binomial generalized linear model with logit link  #
######################################################

library(drc) # to load the data set
data(selenium)
SE <- selenium
SE$t <- factor(SE$type)
SE$typef <- SE$t
levels(SE$typef) <- c("Selenate","Selenite","Selenomethionine","Selenocysteine")
SE$lc10 <- log10(SE$conc)
SE$p <- SE$dead/SE$total
SEn0 <- subset(SE, conc!=0)

fitSE <- glm(cbind(dead, total-dead) ~ typef-1 + lc10:typef,
 data=SEn0, family=binomial())

# Matrices for all pairwise comparisons
CTap  <-  contrMat(n=table(SEn0$typef), type="Tukey")
rownames(CTap)  <-  c("Selenite / Selenate", "Selenomethionine / Selenate",
 "Selenocysteine / Selenate", "Selenomethionine / Selenite", 
 "Selenocysteine / Selenite", "Selenocysteine / Selenomethionine")

# Grid of covariate values
VX <- seq(0.7,2.9,0.2)
ITx <- matrix(rep(1,length(VX)),ncol=1)
XM <- cbind(ITx, VX)
CMap <- kronecker(XM,CTap)
CMap

# Compute and transform the intervals
COMPap <- glht(fitSE, linfct=CMap)
SCIap <- confint(COMPap)
ESCIap <- exp(SCIap$confint)

# Plot corresponding to Figure 8)
dtrtap  <-  factor(rep(rownames(CTap), times=length(VX)), levels=rownames(CTap))
dvxap  <-  rep(10^VX, each=nrow(CTap))
dap <- data.frame(ESCIap, comp=dtrtap, conc=dvxap)

ggplot(dap, aes(x=conc, y=Estimate, ymin=lwr, ymax=upr)) + 
 geom_errorbar(width=0.1) + geom_line() + geom_point(shape=15) + 
 facet_wrap(~comp) + xlab("Concentration") + ylab("Odds \n P(dead)/P(surviving)") +
 geom_hline(yintercept=1) + scale_y_log10() + scale_x_log10()

