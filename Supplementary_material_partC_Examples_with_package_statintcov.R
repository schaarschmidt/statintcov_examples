
# Illustrating the use of package statint (which relies on packages 
# mvtnorm, multcomp, mratios, ggplot2) for reproducing the three examples.
# Packages MASS to load the data of Example 4.1 and package drc to load
# the data of Example 4.3.

library("ggplot2")
library("multcomp")
library("mratios")
library("MASS")
library("drc")

# Installation of package statintcov
# from GitHub, using the package devtools:

install.packages("devtools")
library("devtools")
install_github(repo="schaarschmidt/statintcov")

library("statintcov")

###############################################
# Example 4.1: All pairwise comparisons in    #
# a linear model with baseline as covariate   #
###############################################

# load data from package MASS:

data(anorexia, package="MASS")
AN <- anorexia

fitAN <- lm(Postwt ~ Treat*Prewt, data=AN)
anova(fitAN)

# Simultaneous confidence intervals for the model predictions
# over a grid of covariate Prewt = 75, 80, ..., 90

dpredAN <- expand.grid(Treat=levels(AN$Treat),
                      Prewt=seq(from=70, to=95, by=5))

lnfctpredAN <- model.matrix(~ Treat*Prewt, data=dpredAN)
lnfctpredAN

# The above linear combintaions can be directly assigned
# to argument linfct in glht

scipredAN <- confint(glht(fitAN, linfct=lnfctpredAN))
dscipredAN <- cbind(dpredAN, scipredAN$confint)

# Reproducing Figure 4:

ggplot(data=AN, aes(y=Postwt, x=Prewt)) + geom_point() +  facet_grid(.~Treat) +
  geom_errorbar(data=dscipredAN, aes(y=Estimate, ymin=lwr, ymax=upr), colour="darkgrey", width=1) +
  geom_line(data=dscipredAN, colour="darkgrey", aes(y=Estimate)) + 
  geom_point(data=dscipredAN, aes(y=Estimate), colour="darkgrey", shape=15) +
  xlab("Preweight") + ylab("Expected postweight")


# All pairwise comparisons for a grid of covariate values preweight = 75,76,...,95

dsciAN <- scitreatcov(response="Postwt", treatment="Treat", covariate="Prewt", data=anorexia,
 covset=seq(from=75, to=95, by=5), treatcon="Tukey", conf.level=0.95, base=2)

str(dsciAN, max.level=1)
str(dsciAN$sci)

# Reproducing Figure 5:

 ggplot(dsciAN$sci, aes(x=covariate, y=Estimate, ymin=lwr, ymax=upr)) + 
 geom_errorbar(width=1) + geom_line() + geom_point(shape=15) +
 facet_grid(.~comparison) + xlab("Preweight") + 
 ylab("Difference in expected postweight") + geom_hline(yintercept=0)


###############################################
# Example 4.2: Comparisons to control with    #
# interaction to a quadratic regression term  #
###############################################

# Package statint contains a data set called pc, that resembles the structure of the data set presented
# in Milliken and Johnson (2002). However, the data values are simulated, thus the following code will NOT exactly
# reproduce the estimates and confidence limits shown in the manuscript for this example. 
# Reference: Milliken GA and Johnson DE (2002): Analysis of Messy Data, Vol.II: Analysis of Covariance. Chapman & Hall/HRC.

data(pc)

# Model selected by Milliken and Johnson:

fitPC <- lm(yield ~ x + I(x^2) + additive + additive:I(x^2), data=pc)
anova(fitPC)

# Since this model contains a quadratic term, it can only be handled via the function cmiacov
# 1) computing the contrast matrix, suitable for parameterization in fitPC
# the only covariate in the data is x, although it appears in several terms
# comparsisons to control

cmPC <- cmiacov(fitPC, treatment="additive",
 covset=list("x" = seq(from=0, to=10, by=1)), treatcon="Dunnett")

str(cmPC, max.level=1)
str(cmPC$linfct)

# the matrix of coefficients in cmpc$linfct
# can be passed to the linfct argument of glht (package multcomp)

sciPC <- confint(glht(fitPC, linfct=cmPC$linfct))

# Plot corresponding to Figure 8:
# For plotting: write the comparsions and covariate values in a data.frame
# together with the confidence interavls from glht

dsciPC <- cbind(cmPC$datadiff, sciPC$confint)

ggplot(dsciPC, aes(x=x, y=Estimate, ymin=lwr, ymax=upr)) +
 geom_errorbar(width=1) + geom_line() + geom_point(shape=15) + 
 facet_grid(.~comparison) + xlab("Compound x") + 
 ylab("Difference in expected yields") + geom_hline(yintercept=0)


# the same approach for ratios to control:

cmrPC <- cmratioiacov(fitPC, treatment="additive",
                covset=list("x" = seq(from=1, to=10, by=1)), treatcon="Dunnett")

str(cmrPC, max.level=1)

# the elements numC and denC can be supplied to the function gsci.ratio
# and define the numerators and denominators of the matrices of
# interest


scirPC <- gsci.ratio(est = coef(fitPC), vcmat = vcov(fitPC),
  degfree=fitPC$df.residual, Num.Contrast=cmrPC$numC,
  Den.Contrast=cmrPC$denC, adjusted = TRUE)

scirPC

dscirPC <- cbind(cmrPC$dataratio, scirPC$estimate, scirPC$conf.int)

# Figure 9:

ggplot(dscirPC, aes(x=x, y=estimate, ymin=lower, ymax=upper)) +
 geom_errorbar(width=1) + geom_line() + geom_point(shape=15) + 
 facet_grid(.~comparison) + xlab("Compound x") + 
 ylab("Ratio in expected yields") + geom_hline(yintercept=1)




######################################################
# Example 4.3: All pairwise comparisons in a         #
# binomial generalized linear model with logit link  #
######################################################

data(selenium, package="drc")

SE <- subset(selenium, conc!=0)
SE$typef <- factor(SE$type)
levels(SE$typef) <- c("Selenate","Selenite","Selenomethionine","Selenocysteine")
SE$l10conc <- log10(SE$conc)

SE$mortality <- SE$dead/SE$total

ggplot(SE, aes(y=mortality, x=conc)) + geom_point() + facet_wrap( ~typef) + scale_x_log10()

fitSE <- glm(cbind(dead, total-dead) ~ typef*l10conc, data=SE, family=binomial())
anova(fitSE, test="Chisq")

# the binomial assumption, or logit-log linearity may not be adequat
# using the quasibinomial assumption should be more reliable

# Computing the contrast matrix, suitable for parameterization in fitSe
# all pairwiSe comparisons ("Tukey")
# a Sequence of 7 points over the obServed range of the covariate

range(SE$l10conc)

cmSE <- cmiacov(fitSE, treatment="typef", 
covset=list("l10conc" = seq(from=0.7,to=2.9,by=0.2)), treatcon="Tukey")

str(cmSE, max.level=1)
str(cmSE$linfct)


# Pass the coeffiecients defining the linear combinations in cmSE$linfct
# to the linfct argument of glht

sciSE <- confint(glht(fitSE, linfct=cmSE$linfct))

# Write the comparsions and covariate values in a data.frame
# together with the confidence interavls from glht (for plotting)

dsciSE <- cbind(cmSE$datadiff, sciSE$confint)

# Plot corresponding to Figure 11)

ggplot(dsciSE, aes(x=l10conc, y=exp(Estimate), ymin=exp(lwr), ymax=exp(upr))) +
geom_errorbar(width=0.1) + geom_line() + geom_point(shape=15) +
facet_wrap(~comparison) + xlab("Log10 concentration") +
ylab("Difference (logit-scale)") + geom_hline(yintercept=0) +
 geom_hline(yintercept=1) + scale_y_log10()
# scale_y_log10(breaks=rep(c(1,2,5), times=4)*rep(c(0.01,0.1,1,10), each=3))

# Note that this example data set exhibits substantial ovferdispersion
# and thus analysis under the quasibinomial assumption is much more reliable
