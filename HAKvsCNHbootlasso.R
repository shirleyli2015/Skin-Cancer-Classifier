library(dplyr)
library(magrittr)
library(readxl)
library(ggbiplot)
library(glmnet)
library(pROC)
library(boot)
source('preprocfunctions.R')
source('johnfuncs.R')
data <- read_excel("HAK-v-CNH_objective.xlsx")
yall <- data$Class
yall <- recode(yall, "HAK"='1', "CNH"='0')
xall <- as.matrix(data[,-(1:2)])
booty<-bootLasso(xall, yall, B= 500, nfolds = 5)
pred <- xall[,13]
cm <- table(True=yall, Predict=!pred)
calcAccuracy(cm)
library(xlsx)
write.csv(t(rbind(colnames(data),booty$Beta,booty$interval)), "bootcoefs.csv")

calcAccuracy <- function(cm) {
  classes <- row.names(cm)
  recalls <- sapply(1:length(classes), function(x) cm[x,x]/sum(cm[x,])*100)
  #recall <- data.frame(recall1 = recalls[1], recall2 = recalls[2], accuracy = sum(diag(cm))/sum(cm)*100)
  recall <- c(recall1 = recalls[1], recall2 = recalls[2], accuracy = sum(diag(cm))/sum(cm)*100)
  
  return(recall)
}

statistic = function(data, selection) {
  #Create logistic regression model
  yall <- data$Class[selection]
  xall <- as.matrix(data[selection,-(1:2)])
  fit <- glmnet(xall, yall, family="binomial", alpha =1, lambda.min.ratio = 1e-05)
  cvfit <- cv.glmnet(xall,yall, nfolds= 5, type.measure="class", keep = TRUE, family="binomial", alpha = 1,lambda.min.ratio=1e-05)
  
  minLamIdx <- 
    which(cvfit$lambda == cvfit$lambda.min)
  rocCurve <- roc(yall, cvfit$fit.preval[,minLamIdx])
  optCut <- coords(rocCurve, "best", best.method="youden", ret="threshold",transpose = TRUE)
  #plot(rocCurve, print.thres="best", print.thres.best.method="youden", print.auc=TRUE, auc.polygon=TRUE)
  
  allPred01 <- cvfit$fit.preval > optCut[1]
  classPred <- as.integer(allPred01[, minLamIdx])
  cm <- table(True=yall, Predict=classPred)
  nonzeroIdx <- unlist(predict(fit, s = cvfit$lambda[minLamIdx], type = "nonzero"))
  coef <- unlist(predict(fit, s = cvfit$lambda[minLamIdx], type = "coefficients")[-1])
  #Create ROC curve and return AUC
  return(
    c(calcAccuracy(cm),exp(coef))
  )
}
set.seed(1)
coochieexp <- boot(data, statistic, R=10000)
coochieexp
#indx <- 1:length(coochie$t0)
#cis<- sapply(indx, function(x) boot.ci(coochie, index = x)$normal)
#cma <- vector("list", length = length(coochie$t0))
#cma <- lapply(cma, function(x))
