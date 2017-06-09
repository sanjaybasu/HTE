rm(list=ls())
ptm <- proc.time()
library(matrixStats)
library(psych)
library(bindata)
library(clusterGeneration)
library(MASS)
library(OptimalCutpoints)
library(ResourceSelection)
library(stringr)
library(plyr)
library(copula)
library(glmnet)
library(caret)
library(SuperLearner)
library(cvAUC)
library(ROCR)
library(gplots)
library(randomForest)
library(ggplot2)
library(xgboost)
library(data.table)
library(h2o)  
library(h2oEnsemble)

n = 1e6
p = 20 # total covariates eligible for model selection
preal = 8 # covariates affecting baseline risk
phtepos = 2 # those covariates creating more than average treatment effect (increasing drug benefits)
phteneg = 2 # those covariates creating less than average treatment effect (reducing drug benefits)
peff = phtepos+phteneg
trialpop = 6000
iters=1
case = 1
coriter=1

#for (coriter in 1:4)  {  # iterate the degree of correlation among candidate covariates
   corj=log(coriter)*.45272+.19628  # varies the correlation among covariates stepwise from 0.125 to 0.33 to 0.5 to 0 (if coriter = 4)
  x = round(rCopula(n,normalCopula(corj,dim=p)))  
  if (coriter==4) x = matrix(rbinom(n*p,1,.5),ncol=p)
  #cor(x)

#for (case in 1:4) {  # case (1): no ATE, +/- HTE; (2): +ATE, +/- HTE; (3) +ATE, no HTE; (4) no ATE, no HTE
  
    ate=0*(case==1)+0.3*(case>1)  # chosen to have typical ATE in CVD treatment trials of about a 5% absolute risk reduction
    hte=2*(case<3)+0*(case==3)   # chosen to have typical theorized HTE in CVD treatment trials of about +/- 5% absolute risk IQR
    if (case==4) ate = 0
    if (case==4) hte = 0
    baserisk = logistic(rowSums(x[,1:preal])-7)+rnorm(n,mean = 0, sd = 0.001)  # baseline risk ~2.5% per year x 5 year trial, typical of CVD trials
    baserisk[baserisk<0]=0
    baserisk[baserisk>1]=1
    newrisk = logistic((1-ate)*rowSums(x[,1:preal])-hte*rowSums(x[,1:phtepos])+hte*rowSums(x[,(1+phtepos):peff])-7)+rnorm(n,mean = 0, sd = 0.001)
    newrisk[newrisk<0]=0
    newrisk[newrisk>1]=1
    arr=baserisk-newrisk
    # hist(arr)
    # summary(arr)
    # summary(baserisk)
    # summary(newrisk)
    
    treatment = rbinom(n*p,1,.5)
    y = rbinom(n,1,baserisk*(1-treatment)+newrisk*treatment)
    interactions = as.matrix(x*treatment,ncol=p)
    popdata = data.frame(y,treatment,x,interactions)
    names(popdata) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
    interactionsnorx = as.matrix(x*0,ncol=p)
    popdatanorx = data.frame(y,rep(0,n),x,interactionsnorx)
    names(popdatanorx) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
    interactionsallrx = as.matrix(x*1,ncol=p)
    popdataallrx = data.frame(y,rep(1,n),x,interactionsallrx)
    names(popdatanorx) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
    
    
    #### matrix preallocation ####
    correctharmv = rep(0,iters)
    correctnonev = rep(0,iters)
    correctbenv = rep(0,iters)
    wrongharmv = rep(0,iters)
    wrongnonv = rep(0,iters)
    wrongbenv = rep(0,iters)
    biasv = matrix(0,ncol=trialpop,nrow=iters)
    aucv = rep(0,iters)
    hlpv = rep(0,iters)
    hlchiv = rep(0,iters)
    dummypassv = rep(0, iters)
    dummypassv7 = rep(0, iters)
    dummypassv8 = rep(0, iters)
    
    rfcorrectharmv = rep(0,iters)
    rfcorrectnonev = rep(0,iters)
    rfcorrectbenv = rep(0,iters)
    rfwrongharmv = rep(0,iters)
    rfwrongnonv = rep(0,iters)
    rfwrongbenv = rep(0,iters)
    rfbiasv = matrix(0,ncol=trialpop,nrow=iters)
    rfaucv = rep(0,iters)
    rfhlpv = rep(0,iters)
    rfhlchiv = rep(0,iters)
    rfdummypassv = rep(0, iters)
    rfdummypassv7 = rep(0, iters)
    rfdummypassv8 = rep(0, iters)
    
    valcorrectharmv = rep(0,iters)
    valcorrectnonev = rep(0,iters)
    valcorrectbenv = rep(0,iters)
    valwrongharmv = rep(0,iters)
    valwrongnonv = rep(0,iters)
    valwrongbenv = rep(0,iters)
    valbiasv = matrix(0,ncol=trialpop,nrow=iters)
    valaucv = rep(0,iters)
    valhlpv = rep(0,iters)
    valhlchiv = rep(0,iters)
    valdummypassv = rep(0, iters)
    
    valrfcorrectharmv = rep(0,iters)
    valrfcorrectnonev = rep(0,iters)
    valrfcorrectbenv = rep(0,iters)
    valrfwrongharmv = rep(0,iters)
    valrfwrongnonv = rep(0,iters)
    valrfwrongbenv = rep(0,iters)
    valrfbiasv = matrix(0,ncol=trialpop,nrow=iters)
    valrfaucv = rep(0,iters)
    valrfhlpv = rep(0,iters)
    valrfhlchiv = rep(0,iters)
    valrfdummypassv = rep(0, iters)
    
    
    
    ftrue <- ""
    nextx <- ""
    if (peff>1) {
      for (ii in 1:(preal-1)) {
        nextx <- paste("x",ii, sep="")
        if (ii==1) {name <- nextx}
        if (ii>1) {name <- c(name, nextx)}
        ftrue <- paste(ftrue, nextx, ", ", sep="")
      }
      ftrue <- paste(ftrue, "x", ii+1, sep="")
    } else if (preal==1) {
      ftrue <- "x1"
    }
    for (ii in 1:p) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
    }
    ftrue=c(unlist(strsplit(ftrue,", ")))
    
    ffalse <- ""
    nextx <- ""
    ii=0
    if ((p-preal)>1) {
      for (ii in (preal+1):(p-1)) {
        nextx <- paste("x",ii,sep="")
        if (ii==1) {name <- nextx}
        if (ii>1) {name <- c(name, nextx)}
        ffalse <- paste(ffalse, nextx, ", ", sep="")
      }
      ffalse <- paste(ffalse, "x", p, sep="")
    } else if ((p-preal)==1) {
      ffalse <- paste("x",preal+1)
    }
    for (ii in (preal+1):p) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
    }
    ffalse=c(unlist(strsplit(ffalse,", ")))
  
    
 for (iter in 1:iters) {
      print(coriter)
      print(case)
      print(iter)
      set.seed(iter)
      trialdata = popdata[sample(nrow(popdata), trialpop), ]
      names(trialdata) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
      set.seed(iter)
      trialdatanorx = popdatanorx[sample(nrow(popdatanorx), trialpop), ]
      names(trialdatanorx) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
      set.seed(iter)
      trialdataallrx = popdataallrx[sample(nrow(popdataallrx), trialpop), ]
      names(trialdataallrx) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
      
      set.seed(iter+1)
      valdata = popdata[sample(nrow(popdata), trialpop), ]
      names(valdata) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
      set.seed(iter+1)
      valdatanorx = popdatanorx[sample(nrow(popdatanorx), trialpop), ]
      names(valdatanorx) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
      set.seed(iter+1)
      valdataallrx = popdataallrx[sample(nrow(popdataallrx), trialpop), ]
      names(valdataallrx) <- c("y","treatment",paste("x", 1:p, sep = ""),paste("x", 1:p, "*treatment",sep = ""))
      
      
      set.seed(iter)
      arrsample = sample(arr, trialpop)
      set.seed(iter+1)
      valarrsample = sample(arr, trialpop)
      
      #### classical approach ####
        classicalmodel = (glm(y~.,family=binomial(),data=trialdata))
        classicalmodelaic = stepAIC(classicalmodel,trace=F)


        baseriskest = predict.glm(classicalmodelaic,newdata=trialdatanorx,type="response")
        newriskest = predict.glm(classicalmodelaic,newdata=trialdataallrx,type="response")

        arrest = baseriskest-newriskest
        bias = arrest-arrsample
        correctharm =sum((arrsample<(-0.01))&(arrest<(-0.01)))
        correctnone = sum((arrsample<0.01)&(arrsample>(-0.01))&(arrest<0.01)&(arrest>(-0.01)))
        correctben = sum((arrsample>0.01)&(arrest>0.01))
        wrongharm = sum((arrsample>(-0.01))&(arrest<(-0.01)))
        wrongnon = sum(((arrsample>0.01)|(arrsample<(-0.01)))&(arrest<0.01)&(arrest>(-0.01)))
        wrongben =  sum((arrsample<0.01)&(arrest>0.01))
        print('ready')

        score = predict(classicalmodelaic,newdata=trialdata,type="response")
        cutpointdata = data.frame(score,trialdata$y)
        names(cutpointdata) <- c("score","y")
        optimal.cutpoint.Youden <- optimal.cutpoints(X = score ~ y,methods="Youden",data=cutpointdata,tag.healthy=0)
        auc=optimal.cutpoint.Youden$Youden$Global$measures.acc$AUC
        hlp=hoslem.test(cutpointdata$y, score)$p.value
        hlchi=hoslem.test(cutpointdata$y, score)$statistic


        sumwrongs = (wrongharm+wrongnon+wrongben)/trialpop
        dummypass = (hlp>.05)&(auc[[1]]>.7)&(sumwrongs>.05)
        dummypass7 = (auc[[1]]>.7)&(sumwrongs>.05)
        dummypass8 = (auc[[1]]>.8)&(sumwrongs>.05)

        correctharmv[iter] = correctharm/trialpop
        correctnonev[iter] = correctnone/trialpop
        correctbenv[iter] = correctben/trialpop
        wrongharmv[iter] = wrongharm/trialpop
        wrongnonv[iter] = wrongnon/trialpop
        wrongbenv[iter] = wrongben/trialpop
        biasv[iter,] = bias
        aucv[iter] = auc[1]
        hlpv[iter] = hlp
        hlchiv[iter] = hlchi
        dummypassv[iter] = dummypass
        dummypassv7[iter] = dummypass7
        dummypassv8[iter] = dummypass8


        valbaseriskest = predict.glm(classicalmodelaic,newdata=valdatanorx,type="response")
        valnewriskest = predict.glm(classicalmodelaic,newdata=valdataallrx,type="response")

        valarrest = valbaseriskest-valnewriskest
        valbias = valarrest-valarrsample
        valcorrectharm =sum((valarrsample<(-0.01))&(valarrest<(-0.01)))
        valcorrectnone = sum((valarrsample<0.01)&(valarrsample>(-0.01))&(valarrest<0.01)&(valarrest>(-0.01)))
        valcorrectben = sum((valarrsample>0.01)&(valarrest>0.01))
        valwrongharm = sum((valarrsample>(-0.01))&(valarrest<(-0.01)))
        valwrongnon = sum(((valarrsample>0.01)|(valarrsample<(-0.01)))&(valarrest<0.01)&(valarrest>(-0.01)))
        valwrongben =  sum((valarrsample<0.01)&(valarrest>0.01))

        valscore = predict(classicalmodelaic,newdata=valdata,type="response")
        valcutpointdata = data.frame(valscore,valdata$y)
        names(valcutpointdata) <- c("valscore","y")
        valoptimal.cutpoint.Youden <- optimal.cutpoints(X = valscore ~ y,methods="Youden",data=valcutpointdata,tag.healthy=0)
        valauc=valoptimal.cutpoint.Youden$Youden$Global$measures.acc$AUC
        valhlp=hoslem.test(valcutpointdata$y, valscore)$p.value
        valhlchi=hoslem.test(valcutpointdata$y, valscore)$statistic
        valsumwrongs = (valwrongharm+valwrongnon+valwrongben)/trialpop
        valdummypass = (valhlp>0.05)&(valauc[[1]]>.7)&(valsumwrongs>.05)

        valcorrectharmv[iter] = valcorrectharm/trialpop
        valcorrectnonev[iter] = valcorrectnone/trialpop
        valcorrectbenv[iter] = valcorrectben/trialpop
        valwrongharmv[iter] = valwrongharm/trialpop
        valwrongnonv[iter] = valwrongnon/trialpop
        valwrongbenv[iter] = valwrongben/trialpop
        valbiasv[iter,] = valbias
        valaucv[iter] = valauc[1]
        valhlpv[iter] = valhlp
        valhlchiv[iter] = valhlchi
        valdummypassv[iter] = valdummypass



      #### superlearner ####

     localH2O <- h2o.init(nthreads=-1, max_mem_size="16g") 
        
      train <- as.h2o(trialdata[,1:(p+2)])
      test <- as.h2o(valdata[,1:(p+2)])
      y <- "y"
      x <- setdiff(names(train), y)
      family <- "binomial"
      
      #For binary classification, response should be a factor
      train[,y] <- as.factor(train[,y])  
      test[,y] <- as.factor(test[,y])
      
      
      
      
       learner <- c("h2o.glm.1","h2o.glm.2", "h2o.glm.3", "h2o.randomForest.1", "h2o.gbm.1","h2o.deeplearning.1", "h2o.deeplearning.2")
      #learner <- c("h2o.deeplearning.1", "h2o.deeplearning.2")
      
      metalearner <- "h2o.glm_nn"
      
      family= "binomial"
      
      
      
      
      create_h2o_glm_wrappers <- function(alphas = c(0.0, 0.5, 1.0)) {
        for (i in seq(length(alphas))) {
          alpha <- alphas[i]
          body <- sprintf(' <- function(..., alpha = %s) h2o.glm.wrapper(..., alpha = alpha)', alpha)
          eval(parse(text = paste('h2o.glm.', i, body, sep = '')), envir = .GlobalEnv)
        }
      }
      create_h2o_glm_wrappers()
      h2o.gbm.1 <- function(..., ntrees = 100, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, seed = seed)
      h2o.randomForest.1 <- function(..., ntrees = 100, nbins = 10, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
      h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
      h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
      
      
      h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)
      
      
      rffit = h2o.ensemble(x = 2:(p+2), y = 1, training_frame = train,
                           family = family,
                           learner = learner,
                           metalearner = metalearner,
                           cvControl = list(V = 10, shuffle = TRUE))
      

      
      
      # 
      # rffit = h2o.ensemble(x = 2:(p+2), y = 1, training_frame = train,
      #                      family = family,
      #                      learner = learner,
      #                      metalearner = metalearner,
      #                      cvControl = list(V = 10, shuffle = TRUE))
      # 

      
      #pred <- predict(object = rffit, newdata = train)$pred[,3]
      
      
      # rffit = h2o.deeplearning(y = 1, x = 2:(p+2), training_frame = train)
      
      # rffit <- SuperLearner(Y = trialdata[,"y"], X = trialdata[,2:(p+3-1)],
      #                        SL.library = c("SL.knn","SL.glmnet","SL.randomForest","SL.gam","SL.gbm"),
      #                        verbose = F,
      #                        family = binomial(),
      #                        method = "method.AUC")
      
      #rffit <- gbm(as.formula(paste("y~treatment+", paste(c(ftrue,ffalse), collapse="+"))),distribution = "bernoulli", data=trialdata, n.trees=1000, cv.folds=10, n.cores=2, interaction.depth=3, bag.fraction=0.5, shrinkage=0.1)
      #best.iter <- gbm.perf(rffit,method="cv")
      
      # fitControl <- trainControl(method = "cv", number = 10, repeats = 1, search = "random")
      # rffit <- train(x=trialdata[,2:(p+3-1)],y = as.factor(trialdata[,"y"]),method = "xgbTree", trControl = fitControl)
      
      #rffit <- randomForest(x=trialdata[,2:(p+3-1)],y = as.factor(trialdata[,"y"]),  data=trialdata)
      
      
      trialdatanorxh2o = trialdatanorx[,1:(p+2)]
      trialdatanorxh2o[,1] = factor(trialdatanorxh2o[,1])
      testnorx <- as.h2o(trialdatanorxh2o)
      
      trialdataallrxh2o = trialdataallrx[,1:(p+2)]
      trialdataallrxh2o[,1] = factor(trialdataallrxh2o[,1])
      testallrx <- as.h2o(trialdataallrxh2o)
      
      rfbaseriskest = as.data.frame(predict(rffit,newdata = testnorx)$pred[,3])
      rfnewriskest = as.data.frame(predict(rffit,newdata = testallrx)$pred[,3])
      
      rfarrest = rfbaseriskest-rfnewriskest
      
      rfbias = as.data.frame(rfarrest-arrsample)$p1
      
      rfcorrectharm =sum((arrsample<(-0.01))&(rfarrest<(-0.01)))
      rfcorrectnone = sum((arrsample<0.01)&(arrsample>(-0.01))&(rfarrest<0.01)&(rfarrest>(-0.01)))
      rfcorrectben = sum((arrsample>0.01)&(rfarrest>0.01))
      rfwrongharm = sum((arrsample>(-0.01))&(rfarrest<(-0.01)))
      rfwrongnon = sum(((arrsample>0.01)|(arrsample<(-0.01)))&(rfarrest<0.01)&(rfarrest>(-0.01)))
      rfwrongben =  sum((arrsample<0.01)&(rfarrest>0.01))
      
      rfscore = as.data.frame(predict(rffit,newdata = train)$pred[,3])
      rfauc <- h2o.ensemble_performance(rffit, newdata = train)$ensemble@metrics$AUC
      rfhlp=hoslem.test(trialdata$y, rfscore$p1)$p.value
      rfhlchi=hoslem.test(trialdata$y, rfscore$p1)$statistic
      
      rfsumwrongs = (rfwrongharm+rfwrongnon+rfwrongben)/trialpop
      rfdummypass = (rfhlp>.05)&(rfauc[[1]]>.7)&(rfsumwrongs>.05)
      rfdummypass7 = (rfauc[[1]]>.7)&(rfsumwrongs>.05)
      rfdummypass8 = (rfauc[[1]]>.8)&(rfsumwrongs>.05)
      
      rfcorrectharmv[iter] = rfcorrectharm/trialpop
      rfcorrectnonev[iter] = rfcorrectnone/trialpop
      rfcorrectbenv[iter] = rfcorrectben/trialpop
      rfwrongharmv[iter] = rfwrongharm/trialpop
      rfwrongnonv[iter] = rfwrongnon/trialpop
      rfwrongbenv[iter] = rfwrongben/trialpop
      rfbiasv[iter,] = rfbias
      rfaucv[iter] = rfauc[1]
      rfhlpv[iter] = rfhlp
      rfhlchiv[iter] = rfhlchi
      rfdummypassv[iter] = rfdummypass
      rfdummypassv7[iter] = rfdummypass7
      rfdummypassv8[iter] = rfdummypass8
      
      
      vtrialdatah2o = valdata[,1:(p+2)]
      vtrialdatah2o[,1] = factor(valdata[,1])
      vtrain <- as.h2o(vtrialdatah2o)
      
      vtrialdatanorxh2o = valdatanorx[,1:(p+2)]
      vtrialdatanorxh2o[,1] = factor(valdatanorx[,1])
      vtestnorx <- as.h2o(vtrialdatanorxh2o)
      
      vtrialdataallrxh2o = valdataallrx[,1:(p+2)]
      vtrialdataallrxh2o[,1] = factor(valdataallrx[,1])
      vtestallrx <- as.h2o(vtrialdataallrxh2o)
      
      valrfbaseriskest = as.data.frame(predict(rffit,newdata = vtestnorx)$pred[,3])
      valrfnewriskest = as.data.frame(predict(rffit,newdata = vtestallrx)$pred[,3])
      
      valrfarrest = valrfbaseriskest-valrfnewriskest
      valrfbias = as.data.frame(valrfarrest-valarrsample)$p1
      valrfcorrectharm =sum((valarrsample<(-0.01))&(valrfarrest<(-0.01)))
      valrfcorrectnone = sum((valarrsample<0.01)&(valarrsample>(-0.01))&(valrfarrest<0.01)&(valrfarrest>(-0.01)))
      valrfcorrectben = sum((valarrsample>0.01)&(valrfarrest>0.01))
      valrfwrongharm = sum((valarrsample>(-0.01))&(valrfarrest<(-0.01)))
      valrfwrongnon = sum(((valarrsample>0.01)|(valarrsample<(-0.01)))&(valrfarrest<0.01)&(valrfarrest>(-0.01)))
      valrfwrongben =  sum((valarrsample<0.01)&(valrfarrest>0.01))
      
      valrfscore = as.data.frame(predict(rffit,newdata = vtrain)$pred[,3])
      valrfcutpointdata = data.frame(valrfscore,valdata$y)
      names(valrfcutpointdata) <- c("valrfscore","y")
      valrfoptimal.cutpoint.Youden <- optimal.cutpoints(X = valrfscore ~ y,methods="Youden",data=valrfcutpointdata,tag.healthy=0)
      valrfauc=valrfoptimal.cutpoint.Youden$Youden$Global$measures.acc$AUC
      valrfhlp=hoslem.test(valrfcutpointdata$y, valrfscore$p1)$p.value
      valrfhlchi=hoslem.test(valrfcutpointdata$y, valrfscore$p1)$statistic
      valrfsumwrongs = (valrfwrongharm+valrfwrongnon+valrfwrongben)/trialpop
      valrfdummypass = (valrfhlp>0.05)&(valrfauc[[1]]>.7)&(valrfsumwrongs>.05)
      
      valrfcorrectharmv[iter] = valrfcorrectharm/trialpop
      valrfcorrectnonev[iter] = valrfcorrectnone/trialpop
      valrfcorrectbenv[iter] = valrfcorrectben/trialpop
      valrfwrongharmv[iter] = valrfwrongharm/trialpop
      valrfwrongnonv[iter] = valrfwrongnon/trialpop
      valrfwrongbenv[iter] = valrfwrongben/trialpop
      valrfbiasv[iter,] = valrfbias
      valrfaucv[iter] = valrfauc[1]
      valrfhlpv[iter] = valrfhlp
      valrfhlchiv[iter] = valrfhlchi
      valrfdummypassv[iter] = valrfdummypass
      
      
      h2o.shutdown(prompt = F)
   }
    
    
    statsout = matrix(c(mean(biasv),quantile(biasv,c(.025,.975)),
                        mean(rfbiasv),quantile(rfbiasv,c(.025,.975)),
                        mean(aucv),quantile(aucv,c(.025,.975)),
                        mean(rfaucv),quantile(rfaucv,c(.025,.975)),
                        mean(hlpv),quantile(hlpv,c(.025,.975)),
                        mean(rfhlpv),quantile(rfhlpv,c(.025,.975)),
                        mean(valbiasv),quantile(valbiasv,c(.025,.975)),
                        mean(valrfbiasv),quantile(valrfbiasv,c(.025,.975)),
                        mean(valaucv),quantile(valaucv,c(.025,.975)),
                        mean(valrfaucv),quantile(valrfaucv,c(.025,.975)),
                        mean(valhlpv),quantile(valhlpv,c(.025,.975)),
                        mean(valrfhlpv),quantile(valrfhlpv,c(.025,.975))),ncol=3,byrow=T)
    
    clinout=matrix(c(mean(correctharmv),quantile(correctharmv,c(.025,.975)),mean(rfcorrectharmv),quantile(rfcorrectharmv,c(.025,.975)),
                     mean(correctnonev),quantile(correctnonev,c(.025,.975)),mean(rfcorrectnonev),quantile(rfcorrectnonev,c(.025,.975)),
                     mean(correctbenv),quantile(correctbenv,c(.025,.975)),mean(rfcorrectbenv),quantile(rfcorrectbenv,c(.025,.975)),
                     mean(wrongharmv),quantile(wrongharmv,c(.025,.975)), mean(rfwrongharmv),quantile(rfwrongharmv,c(.025,.975)),
                     mean(wrongnonv),quantile(wrongnonv,c(.025,.975)),mean(rfwrongnonv),quantile(rfwrongnonv,c(.025,.975)),
                     mean(wrongbenv),quantile(wrongbenv,c(.025,.975)),mean(rfwrongbenv),quantile(rfwrongbenv,c(.025,.975)),
                     mean(valcorrectharmv),quantile(valcorrectharmv,c(.025,.975)),mean(valrfcorrectharmv),quantile(valrfcorrectharmv,c(.025,.975)),
                     mean(valcorrectnonev),quantile(valcorrectnonev,c(.025,.975)),mean(valrfcorrectnonev),quantile(valrfcorrectnonev,c(.025,.975)),
                     mean(valcorrectbenv),quantile(valcorrectbenv,c(.025,.975)),mean(valrfcorrectbenv),quantile(valrfcorrectbenv,c(.025,.975)),
                     mean(valwrongharmv),quantile(valwrongharmv,c(.025,.975)),mean(valrfwrongharmv),quantile(valrfwrongharmv,c(.025,.975)),
                     mean(valwrongnonv),quantile(valwrongnonv,c(.025,.975)),mean(valrfwrongnonv),quantile(valrfwrongnonv,c(.025,.975)),
                     mean(valwrongbenv),quantile(valwrongbenv,c(.025,.975)),mean(valrfwrongbenv),quantile(valrfwrongbenv,c(.025,.975))),ncol=6,byrow=T)
    
    errcon = colSums(matrix(c(dummypassv,dummypassv7,dummypassv8,rfdummypassv,rfdummypassv7,rfdummypassv8,valdummypassv,valrfdummypassv),ncol=8,byrow=T))/iters
    
    statsout = data.frame(statsout,row.names=c("bias conv", "bias gbm","c stat conv","c stat gbm",  "gnd p conv","gnd p gbm","val bias conv", "val bias gbm", "val c stat conv", "val c stat gbm","val gnd p conv","val gnd p gbm"))
    colnames(statsout)=c("mean"," 95% low"," 95% high")
    
    clinout = data.frame(clinout,row.names=c("correct harm prediction", "correct neut prediction","correct ben prediction","wrong harm prediction", "wrong neut prediction","wrong ben prediction","val correct harm prediction","val correct neut prediction", "val correct ben prediction","val wrong harm prediction","val wrong neut prediction", "val wrong ben prediction"))
    colnames(clinout)=c("conventional"," 95% low"," 95% high", "gbm"," 95% low"," 95% high")
    
    errcon = data.frame(errcon,row.names=c("conventional calibration","conventional 7","conventional 8","gbm calibration","gbm 7","gbm 8","conventional validation","gbm validation"))
    colnames(errcon)=c(">5% incorrectly predicted, Cstat>0.7, HLtest pass")

    save(statsout, file=paste("HTEstatsout",case,coriter,".RData",sep=""))
    save(clinout, file=paste("HTEclinout",case,coriter,".RData",sep=""))
    save(errcon, file=paste("HTEerrcon",case,coriter,".RData",sep=""))


#  }
  
#}

proc.time() - ptm


