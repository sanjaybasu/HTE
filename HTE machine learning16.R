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
iters=1000

for (coriter in 1:4)  {  # iterate the degree of correlation among candidate covariates
   corj=log(coriter)*.45272+.19628  # varies the correlation among covariates stepwise from 0.125 to 0.33 to 0.5 to 0 (if coriter = 4)
  x = round(rCopula(n,normalCopula(corj,dim=p)))  
  if (coriter==4) x = matrix(rbinom(n*p,1,.5),ncol=p)
  #cor(x)



for (case in 1:4) {  # case (1): no ATE, +/- HTE; (2): +ATE, +/- HTE; (3) +ATE, no HTE; (4) no ATE, no HTE
  
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
    cordummypassv = rep(0,iters)
    sumwrongsv = rep(0,iters)
    corgoodv = rep(0,iters)
    loosecordummypassv = rep(0,iters)
    loosecorgoodv = rep(0,iters)
    
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
    rfcordummypassv = rep(0,iters)
    rfsumwrongsv  = rep(0,iters)
    rfcorgoodv = rep(0,iters)
    rfloosecordummypassv = rep(0,iters)
    rfloosecorgoodv = rep(0,iters)
    
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
        
        newrisktrial = newrisk[sample(nrow(popdata), trialpop)]
        baserisktrial = baserisk[sample(nrow(popdata), trialpop)]
        
        
        bencats = c(-.01,.01)
        
        bencat = 1*(arrest<=bencats[1])+
          2*((arrest>bencats[1])&(arrest<bencats[2]))+
          3*((arrest>=bencats[2]))
        
        correstab = describeBy(trialdata$y,list(bencat,trialdata$treatment),mat=TRUE)
        test1=prop.test(x=c(correstab[1,5]*correstab[1,6],correstab[4,5]*correstab[4,6]), n=c(correstab[1,5],correstab[4,5]), correct=T)
        test2=prop.test(x=c(correstab[2,5]*correstab[2,6],correstab[5,5]*correstab[5,6]), n=c(correstab[2,5],correstab[5,5]), correct=T)
        test3=prop.test(x=c(correstab[3,5]*correstab[3,6],correstab[6,5]*correstab[6,6]), n=c(correstab[3,5],correstab[6,5]), correct=T)
      
        correstest1 = (test1$estimate["prop 2"]<test1$estimate["prop 1"])&(test1$p.value<.05)
        correstest2 = (test2$p.value>0.05)
        correstest3 = (test3$estimate["prop 2"]>test3$estimate["prop 1"])&(test3$p.value<.05)
        correstest = (correstest1==1)&(correstest2==1)&(correstest3==1)
        
        cordummypass = (sumwrongs>0.05)&(correstest==1)
        corgood = (sumwrongs<0.05)&(correstest==1)

        correstest1 = (test1$estimate["prop 2"]<test1$estimate["prop 1"])
        correstest2 = (test2$p.value>0.05)
        correstest3 = (test3$estimate["prop 2"]>test3$estimate["prop 1"])
        correstest = (correstest1==1)&(correstest2==1)&(correstest3==1)
        
        loosecordummypass = (sumwrongs>0.05)&(correstest==1)
        loosecorgood = (sumwrongs<0.05)&(correstest==1)
        

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
        cordummypassv[iter] = cordummypass
        sumwrongsv[iter] = sumwrongs
        corgoodv[iter] = corgood
        loosecordummypassv[iter] = loosecordummypass
        loosecorgoodv[iter] = loosecordummypasscorgood
        
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
        h2o.removeAll() 
        
      train <- as.h2o(trialdata[,1:(p+2)])
      test <- as.h2o(valdata[,1:(p+2)])
      y <- "y"
      x <- setdiff(names(train), y)
      family <- "binomial"
      
      #For binary classification, response should be a factor
      train[,y] <- as.factor(train[,y])  
      test[,y] <- as.factor(test[,y])
      
      
      
      # Random Grid Search (e.g. 120 second maximum)
      # This is set to run fairly quickly, increase max_runtime_secs 
      # or max_models to cover more of the hyperparameter space.
      # Also, you can expand the hyperparameter space of each of the 
      # algorithms by modifying the hyper param code below.
      
      search_criteria <- list(strategy = "RandomDiscrete", 
                              max_runtime_secs = 120)
      nfolds <- 10
      
      # GBM Hyperparameters
      learn_rate_opt <- c(0.01, 0.03) 
      max_depth_opt <- c(3, 4, 5, 6, 9)
      sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
      col_sample_rate_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
      hyper_params <- list(learn_rate = learn_rate_opt,
                           max_depth = max_depth_opt, 
                           sample_rate = sample_rate_opt,
                           col_sample_rate = col_sample_rate_opt)
      
      gbm_grid <- h2o.grid("gbm", x = x, y = y,
                           training_frame = train,
                           ntrees = 2000,
                           seed = 1,
                           nfolds = nfolds,
                           fold_assignment = "Modulo",
                           keep_cross_validation_predictions = TRUE,
                           hyper_params = hyper_params,
                           search_criteria = search_criteria)
      gbm_models <- lapply(gbm_grid@model_ids, function(model_id) h2o.getModel(model_id))
      
      
      
      # RF Hyperparamters
      mtries_opt <- 8:20 
      max_depth_opt <- c(5, 10, 15, 20, 25)
      sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
      col_sample_rate_per_tree_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
      hyper_params <- list(mtries = mtries_opt,
                           max_depth = max_depth_opt,
                           sample_rate = sample_rate_opt,
                           col_sample_rate_per_tree = col_sample_rate_per_tree_opt)
      
      rf_grid <- h2o.grid("randomForest", x = x, y = y,
                          training_frame = train,
                          ntrees = 2000,
                          seed = 1,
                          nfolds = nfolds,
                          fold_assignment = "Modulo",
                          keep_cross_validation_predictions = TRUE,                    
                          hyper_params = hyper_params,
                          search_criteria = search_criteria)
      rf_models <- lapply(rf_grid@model_ids, function(model_id) h2o.getModel(model_id))
      
      
      
      # Deeplearning Hyperparamters
      activation_opt <- c("Rectifier", "RectifierWithDropout", 
                          "Maxout", "MaxoutWithDropout") 
      hidden_opt <- list(c(10,10), c(20,15), c(50,50,50))
      l1_opt <- c(0, 1e-3, 1e-5)
      l2_opt <- c(0, 1e-3, 1e-5)
      hyper_params <- list(activation = activation_opt,
                           hidden = hidden_opt,
                           l1 = l1_opt,
                           l2 = l2_opt)
      
      dl_grid <- h2o.grid("deeplearning", x = x, y = y,
                          training_frame = train,
                          epochs = 15,
                          seed = 1,
                          nfolds = nfolds,
                          fold_assignment = "Modulo",
                          keep_cross_validation_predictions = TRUE,                    
                          hyper_params = hyper_params,
                          search_criteria = search_criteria)
      dl_models <- lapply(dl_grid@model_ids, function(model_id) h2o.getModel(model_id))
      
      
      
      # GLM Hyperparamters
      alpha_opt <- seq(0,1,0.1)
      lambda_opt <- c(0,1e-7,1e-5,1e-3,1e-1)
      hyper_params <- list(alpha = alpha_opt,
                           lambda = lambda_opt)
      
      glm_grid <- h2o.grid("glm", x = x, y = y,
                           training_frame = train,
                           family = "binomial",
                           nfolds = nfolds,
                           fold_assignment = "Modulo",
                           keep_cross_validation_predictions = TRUE,                    
                           hyper_params = hyper_params,
                           search_criteria = search_criteria)
      glm_models <- lapply(glm_grid@model_ids, function(model_id) h2o.getModel(model_id))
      
      
      # Create a list of all the base models
      models <- c(gbm_models, rf_models, dl_models, glm_models)
      
      # Specify a nonnegative weight GLM as the metalearner
      h2o.glm_nn <- function(..., non_negative = T) h2o.glm.wrapper(..., non_negative = non_negative) # define meta-learner [GLM restricted to non-neg weights, which is shown in the literature to improve outcomes from ensembles]
      
      metalearner <- "h2o.glm_nn"
      
      # Stack models
      stack <- h2o.stack(models = models, 
                         response_frame = train[,y],
                         metalearner = metalearner)
      
      
      # GLM restricted to non-negative weights as a metalearner
      h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)
      rffit <- h2o.metalearn(stack, metalearner = "h2o.glm_nn")
      perf <- h2o.ensemble_performance(rffit, newdata = train, score_base_models = T)
      

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
      
      rfbencat = 1*(rfarrest<=bencats[1])+
        2*((rfarrest>bencats[1])&(rfarrest<bencats[2]))+
        3*((rfarrest>=bencats[2]))
      
      rfcorrestab = describeBy(trialdata$y,list(rfbencat,trialdata$treatment),mat=TRUE)
      rftest1=prop.test(x=c(rfcorrestab[1,5]*rfcorrestab[1,6],rfcorrestab[4,5]*rfcorrestab[4,6]), n=c(rfcorrestab[1,5],rfcorrestab[4,5]), correct=T)
      rftest2=prop.test(x=c(rfcorrestab[2,5]*rfcorrestab[2,6],rfcorrestab[5,5]*rfcorrestab[5,6]), n=c(rfcorrestab[2,5],rfcorrestab[5,5]), correct=T)
      rftest3=prop.test(x=c(rfcorrestab[3,5]*rfcorrestab[3,6],rfcorrestab[6,5]*rfcorrestab[6,6]), n=c(rfcorrestab[3,5],rfcorrestab[6,5]), correct=T)
      
      rfcorrestest1 = (rftest1$estimate["prop 2"]<rftest1$estimate["prop 1"])&(rftest1$p.value<.05)
      rfcorrestest2 = (rftest2$p.value>0.05)
      rfcorrestest3 = (rftest3$estimate["prop 2"]>rftest3$estimate["prop 1"])&(rftest3$p.value<.05)
      rfcorrestest = (rfcorrestest1==1)&(rfcorrestest2==1)&(rfcorrestest3==1)
      
      rfcordummypass = (rfsumwrongs>0.05)&(rfcorrestest==1)
      rfcorgood = (rfsumwrongs<0.05)&(rfcorrestest==1)
      
      rfcorrestest1 = (rftest1$estimate["prop 2"]<rftest1$estimate["prop 1"])
      rfcorrestest2 = (rftest2$p.value>0.05)
      rfcorrestest3 = (rftest3$estimate["prop 2"]>rftest3$estimate["prop 1"])
      rfcorrestest = (rfcorrestest1==1)&(rfcorrestest2==1)&(rfcorrestest3==1)
      
      rfloosecordummypass = (rfsumwrongs>0.05)&(rfcorrestest==1)
      rfloosecorgood = (rfsumwrongs<0.05)&(rfcorrestest==1)
      
      
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
      rfcordummypassv[iter] = rfcordummypass
      rfsumwrongsv[iter] = rfsumwrongs
      rfcorgoodv[iter] = rfcorgood
      rfloosecordummypassv[iter] = rfloosecordummypass
      rfloosecorgoodv[iter] = rfloosecorgood
      
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
      
      
    }
    
    
    statsout = matrix(c(mean(na.omit(biasv)),quantile(biasv,c(.025,.975)),
                        mean(na.omit(rfbiasv)),quantile(rfbiasv,c(.025,.975)),
                        mean(na.omit(aucv)),quantile(aucv,c(.025,.975)),
                        mean(na.omit(rfaucv)),quantile(rfaucv,c(.025,.975)),
                        mean(na.omit(hlpv)),quantile(hlpv,c(.025,.975)),
                        mean(na.omit(rfhlpv)),quantile(rfhlpv,c(.025,.975)),
                        mean(na.omit(valbiasv)),quantile(valbiasv,c(.025,.975)),
                        mean(na.omit(valrfbiasv)),quantile(valrfbiasv,c(.025,.975)),
                        mean(na.omit(valaucv)),quantile(valaucv,c(.025,.975)),
                        mean(na.omit(valrfaucv)),quantile(valrfaucv,c(.025,.975)),
                        mean(na.omit(valhlpv)),quantile(valhlpv,c(.025,.975)),
                        mean(na.omit(valrfhlpv)),quantile(valrfhlpv,c(.025,.975))),ncol=3,byrow=T)
    
    clinout=matrix(c(mean(na.omit(correctharmv)),quantile(correctharmv,c(.025,.975)),mean(na.omit(rfcorrectharmv)),quantile(rfcorrectharmv,c(.025,.975)),
                     mean(na.omit(correctnonev)),quantile(correctnonev,c(.025,.975)),mean(na.omit(rfcorrectnonev)),quantile(rfcorrectnonev,c(.025,.975)),
                     mean(na.omit(correctbenv)),quantile(correctbenv,c(.025,.975)),mean(na.omit(rfcorrectbenv)),quantile(rfcorrectbenv,c(.025,.975)),
                     mean(na.omit(wrongharmv)),quantile(wrongharmv,c(.025,.975)), mean(na.omit(rfwrongharmv)),quantile(rfwrongharmv,c(.025,.975)),
                     mean(na.omit(wrongnonv)),quantile(wrongnonv,c(.025,.975)),mean(na.omit(rfwrongnonv)),quantile(rfwrongnonv,c(.025,.975)),
                     mean(na.omit(wrongbenv)),quantile(wrongbenv,c(.025,.975)),mean(na.omit(rfwrongbenv)),quantile(rfwrongbenv,c(.025,.975)),
                     mean(na.omit(valcorrectharmv)),quantile(valcorrectharmv,c(.025,.975)),mean(na.omit(valrfcorrectharmv)),quantile(valrfcorrectharmv,c(.025,.975)),
                     mean(na.omit(valcorrectnonev)),quantile(valcorrectnonev,c(.025,.975)),mean(na.omit(valrfcorrectnonev)),quantile(valrfcorrectnonev,c(.025,.975)),
                     mean(na.omit(valcorrectbenv)),quantile(valcorrectbenv,c(.025,.975)),mean(na.omit(valrfcorrectbenv)),quantile(valrfcorrectbenv,c(.025,.975)),
                     mean(na.omit(valwrongharmv)),quantile(valwrongharmv,c(.025,.975)),mean(na.omit(valrfwrongharmv)),quantile(valrfwrongharmv,c(.025,.975)),
                     mean(na.omit(valwrongnonv)),quantile(valwrongnonv,c(.025,.975)),mean(na.omit(valrfwrongnonv)),quantile(valrfwrongnonv,c(.025,.975)),
                     mean(na.omit(valwrongbenv)),quantile(valwrongbenv,c(.025,.975)),mean(na.omit(valrfwrongbenv)),quantile(valrfwrongbenv,c(.025,.975))),ncol=6,byrow=T)
    
    newvec = matrix(c(mean(na.omit(sumwrongsv)), quantile(sumwrongsv,c(.025,.975)),
                      mean(na.omit(cordummypassv)), quantile(cordummypassv,c(.025,.975)),
                      mean(na.omit(corgoodv)), quantile(corgoodv, c(.025,.975)),
                      mean(na.omit(loosecordummypassv)), quantile(loosecordummypassv,c(.025,.975)),
                      mean(na.omit(loosecorgoodv)), quantile(loosecorgoodv, c(.025,.975)),
                      mean(na.omit(rfsumwrongsv)), quantile(rfsumwrongsv,c(.025,.975)),
                      mean(na.omit(rfcordummypassv)), quantile(rfcordummypassv,c(.025,.975)),
                      mean(na.omit(rfcorgoodv)), quantile(rfcorgoodv, c(.025,.975)),
                      mean(na.omit(rfloosecordummypassv)), quantile(rfloosecordummypassv,c(.025,.975)),
                      mean(na.omit(rfloosecordummypassv)), quantile(rfloosecordummypassv, c(.025,.975))), ncol=3,byrow=T)
    
    errcon = colSums(matrix(c(dummypassv,dummypassv7,dummypassv8,rfdummypassv,rfdummypassv7,rfdummypassv8,valdummypassv,valrfdummypassv),ncol=8,byrow=T))/iters
    
    statsout = data.frame(statsout,row.names=c("bias conv", "bias gbm","c stat conv","c stat gbm",  "gnd p conv","gnd p gbm","val bias conv", "val bias gbm", "val c stat conv", "val c stat gbm","val gnd p conv","val gnd p gbm"))
    colnames(statsout)=c("mean"," 95% low"," 95% high")
    
    clinout = data.frame(clinout,row.names=c("correct harm prediction", "correct neut prediction","correct ben prediction","wrong harm prediction", "wrong neut prediction","wrong ben prediction","val correct harm prediction","val correct neut prediction", "val correct ben prediction","val wrong harm prediction","val wrong neut prediction", "val wrong ben prediction"))
    colnames(clinout)=c("conventional"," 95% low"," 95% high", "gbm"," 95% low"," 95% high")
    
    newvec = data.frame(newvec,row.names=c("wrong perc conv", "bad corres pass conv","good corres pass conv","bad loose corres pass conv","good loose corres pass conv", "wrong perc ML", "bad corres pass ML", "good corres pass ML", "bad loose corres pass ML", "good loose corres pass ML"))
    colnames(newvec)=c("mean","95lo","95hi")
    
    errcon = data.frame(errcon,row.names=c("conventional calibration","conventional 7","conventional 8","gbm calibration","gbm 7","gbm 8","conventional validation","gbm validation"))
    colnames(errcon)=c(">5% incorrectly predicted, Cstat>0.7, HLtest pass")

    save(statsout, file=paste("HTEstatsout",case,coriter,".RData",sep=""))
    save(clinout, file=paste("HTEclinout",case,coriter,".RData",sep=""))
    save(errcon, file=paste("HTEerrcon",case,coriter,".RData",sep=""))
    save(newvec, file=paste("HTEnewout",case,coriter,".RData",sep=""))


 }

}

proc.time() - ptm


