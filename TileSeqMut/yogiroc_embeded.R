# Title     : TODO
# Objective : TODO
# Created by: roujia
# Created on: 2021-06-29

#load library
library(yogiroc)

args<-commandArgs(TRUE)

inputfile <- args[1]
plotTitle <- args[2]
plotFile <- args[3]
inputData <- read.csv(file = inputfile, colClasses=c("pathogenic"="logical"))

##load data from file
#fakeData <- data.frame(
#  pathogenic = c(rep(TRUE,90),rep(FALSE,60)),
#  prediction1 = c(rnorm(90,mean=2,sd=1),rnorm(60,mean=1,sd=1)),
#  prediction2 = c(rnorm(90,mean=1.5,sd=1),rnorm(60,mean=1,sd=1))
#)
#

#print(inputData[ , -1, drop=FALSE])
#create yogiroc2 object
#print(inputData)
yrobj <- yr2(truth=inputData$pathogenic, scores=inputData[ , -1, drop=FALSE])
#draw PRC curves
#draw.prc(yrobj)

#draw prior-balanced PRC curves
pdf(plotFile)
draw.prc.CI(yrobj,balanced=TRUE,main=plotTitle)
dev.off()
#use custom colors and other graphical parameters
#draw.prc(yrobj,col=c("chartreuse3","firebrick3"),main="My PRC curve")

#draw PRC curves without monotonization
#draw.prc(yrobj,monotonized=FALSE)

#draw PRC curves with confidence itervals
#draw.prc.CI(yrobj,balanced=TRUE,main=plotTitle)

#calculate AUPRC / AUBPRC / R90P
#auprc(yrobj)
#auprc(yrobj,balanced=TRUE)
#recall.at.prec(yrobj,0.9,balanced=TRUE)
#
##calculate whether differences in AUPRC are significant
#aucSignif <- auprc.signif(yrobj)
#sprintf("Prediction 1 has an AUPRC of %.02f with confidence interval [%.02f;%.02f].",
#  aucSignif$auprc[[1]], aucSignif$ci[1,1], aucSignif$ci[2,1]
#)
#sprintf("Prediction 2 has an AUPRC of %.02f with confidence interval [%.02f;%.02f].",
#  aucSignif$auprc[[2]], aucSignif$ci[1,2], aucSignif$ci[2,2]
#)
#sprintf("Prediction 1 is better than Prediction 2, with p-value=%.03f and LLR=%.02f",
#   aucSignif$pval[1,2], aucSignif$llr[1,2]
#)
