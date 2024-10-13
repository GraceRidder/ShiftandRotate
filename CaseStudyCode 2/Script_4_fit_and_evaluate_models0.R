sl = 10 #current slot
nslots = 10

library(Hmsc)

# make a folder in your directory called "Models" 
wd = "~/Desktop/CaseStudyCode2"
model.directory = file.path(wd, "Models")

##### within the thin loop with the rest of the models ######

load("unfitted_models.RData")#models
nm = length(models)

re = rep(NA,nm)
for(i in 1:nm) re[i] = (models[[i]]$nr==0)
for(i in 1:nm) re[i] = TRUE
sel = which(re)
nm = length(sel)

xx = round(seq(from=1,to=nm,length=nslots+1))
starts = xx[1:(nslots)]
ends = xx[2:(nslots+1)]-1
ends[nslots] = xx[nslots+1]


nChains = 2
verbose = 0
thin.list = c(1,1,10,100)
samples.list = c(5,250,250,250)

for(ts in 1:2){
  thin=thin.list[ts]
  samples=samples.list[ts]
  transient = round(samples/2)*thin
  for(i in sel[starts[sl]:ends[sl]]){
    filename=file.path(model.directory, 
                       paste0("results_",as.character(ts),"/",
                              "model_",as.character(i),
                              "_samples_",as.character(samples),
                              "_thin_",as.character(thin),".RData"))
    if(file.exists(filename)){
#      print("file exists:")
      load(filename)
      if(any(is.na(AUC.tot))){
        print(filename)
        print("NA found")
        file.rename(filename,paste0(filename,".old"))
      }
    } else
    {
      print(filename)      
      m = sampleMcmc(models[[i]], thin = thin, samples = samples, 
                     transient = transient, nChains = nChains, verbose = verbose)
      
      ##Convergence statistics 
      mpost = convertToCodaObject(m)
      psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
      BetaGelman <- psrf.beta 
      
      ##inference
      postBeta = getPostEstimate(m, parName="Beta")
      
      # (1) based on fitting the models to the full data (in which case we more precisely measure
      # explanatory rather than predictive power), 
      preds = computePredictedValues(m)
      MF = evaluateModelFit(hM=m, predY=preds) 
      
      #(2) by applying two-fold cross-validation where the sampling units were fully randomly assigned to
      # the two folds (measuring predictive power for spatial and environmental interpolation),
      partition = createPartition(m, nfolds = 2)
      preds = computePredictedValues(m,partition=partition)
      MFCV = evaluateModelFit(hM=m, predY=preds)
      
      #(3) by applying two-fold cross-validation where the sampling units were assigned to the two 
      # folds deterministically based on the longitude (measuring predictive power for spatial extrapolation),
      lat <- rank(m[["XData"]][[2]])
      for (h in 1:length(lat)){
        if (lat[h] > round(length(lat)/2)){
          lat[h] = 2
        } else {
          lat[h] = 1
        }
      }
      
      partition <- as.integer(lat)
      preds = computePredictedValues(m,partition=partition)
      MFCVS = evaluateModelFit(hM=m, predY=preds)
      
      # (4) by applying two-fold cross-validation where the sampling units were assigned to the two folds 
      # deterministically based on the first environmental predictor (measuring predictive power for 
      # environmental extrapolation)
      elv <- rank(m[["XData"]][[3]])
      for (h in 1:length(elv)){
        if (elv[h] > round(length(lat)/2)){
          elv[h] = 2
        } else {
          elv[h] = 1
        }
      }
      
      partition <- as.integer(elv)
      preds = computePredictedValues(m,partition=partition)
      MFCVE = evaluateModelFit(hM=m, predY=preds)
      
      tmp = c(mean(MF$AUC,na.rm = TRUE),
              mean(MFCV$AUC,na.rm = TRUE),
              mean(MFCVS$AUC,na.rm = TRUE),
              mean(MFCVE$AUC,na.rm = TRUE))
      AUC.tot <- tmp
      
      save(BetaGelman,postBeta,AUC.tot,file=filename)
    }
  }
}
