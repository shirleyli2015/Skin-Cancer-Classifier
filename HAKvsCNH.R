library(dplyr)
library(magrittr)
library(readxl)
library(ggbiplot)
library(glmnet)
library(pROC)
source('preprocfunctions.R')
source('johnfuncs.R')
set.seed(5)
data <- read_excel("HAK-v-CNH_objective.xlsx")

yall <- data$Class
xall <- as.matrix(data[,-(1:2)])
xall <- scale(xall, TRUE, TRUE)
yall <- as.numeric(recode(yall, "HAK"='1', "CNH"='0'))
fit <- glmnet(xall, yall, family="binomial", alpha =1, lambda.min.ratio = 1e-05, exact = TRUE, standardize = FALSE)
cvfit <- cv.glmnet(xall,yall, nfolds= 5, type.measure="class", keep = TRUE, family="binomial", alpha = 1,lambda.min.ratio=1e-05)

minLamIdx <- which(cvfit$lambda == cvfit$lambda[24])
  #which(cvfit$lambda == cvfit$lambda.1se)

rocCurve <- roc(yall, cvfit$fit.preval[,minLamIdx])
optCut <- coords(rocCurve, "best", best.method="youden", ret="threshold",transpose = TRUE)
plot(rocCurve, print.thres="best", print.thres.best.method="youden", print.auc=TRUE, auc.polygon=TRUE)

allPred01 <- cvfit$fit.preval > optCut[1]
classPred <- as.integer(allPred01[, minLamIdx])
cm <- table(True=yall, Predict=classPred)
cat("\n")

cm
bob<-calcAccuracy(cm)

optCut
nonzeroIdx <- unlist(predict(fit, s = cvfit$lambda[minLamIdx], type = "nonzero"))
coef <- unlist(predict(fit, s = cvfit$lambda[minLamIdx], type = "coefficients")[-1])
coef <- coef[coef != 0]
sortIdx <- order(coef)
nzMZcoef <- cbind(colnames(xall)[nonzeroIdx], coef)[sortIdx,]
nzMZcoef

plot(cvfit$glmnet.fit, "norm",   label=TRUE)
plot(cvfit$glmnet.fit, "lambda", label=TRUE)
plot(cvfit)

classes <- row.names(cm)
recalls <- sapply(1:length(classes), function(x) cm[x,x]/sum(cm[x,])*100)
cat(sprintf("%s Recall: %.4f%%\n", classes, recalls))
cat(sprintf("Overall Accuracy: %.4f%%\n", sum(diag(cm))/sum(cm)*100))

#fixed Lasso Inf
library(selectiveInference)
beta = coef(fit, s=cvfit$lambda[minLamIdx])[-1]

fixedLassoInf(xall, yall,beta = beta, lambda = cvfit$lambda[minLamIdx]* length(yall))

