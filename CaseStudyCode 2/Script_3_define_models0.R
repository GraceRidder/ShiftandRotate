load("datasets.RData")


library(Hmsc)
set.seed(1)
nd = 100
p = sample(nd)
while(any(p==1:nd)) p = sample(nd)

models = list()
modelnames = c()
ind = 0
for(i in 1:nd){
  X = environmental.covariates[[i]]
  X.among = environmental.covariates[[p[i]]]
  X.within = X[sample(nrow(X)),]
  xy = spatial.coordinates[[i]]
  Ys = data.sets[[i]]
  studyDesign <- data.frame(route=rownames(xy))
  studyDesign$route <-as.factor(studyDesign$route)
  rL = HmscRandomLevel(sData=xy)
  
  for(j in 1:4){
    Y = Ys[[j]]
    for(k in 1:3){
      if(k==1){
        XData = as.data.frame(X)
        dname = names(Ys)[j]
      }
      if(k==2){
        XData = as.data.frame(X.within)
        dname = paste0(names(Ys)[j],".W")
      }
      if(k==3){
        XData = as.data.frame(X.among)
        dname = paste0(names(Ys)[j],".A")
      }
      for(l1 in 1:2){
        if(l1==1){
          XFormula = ~elv+ BIO4 + BIO5 + BIO6 + BIO7 + BIO12
          mname1="M.E6"
        }
        if(l1==2){
          XFormula =  ~elv
          mname1="M.E1"
        }
        for(l2 in 1:2){
          if(l2==1){
            m = Hmsc(Y=Y, XData = XData, XFormula= XFormula,
                     studyDesign = studyDesign, distr = "probit")
            mname = mname1
          }
          if(l2==2){
            m = Hmsc(Y=Y, XData = XData, XFormula= XFormula,
                     studyDesign = studyDesign, ranLevels=list('route'=rL), distr = "probit")
            mname=paste0(mname1,".S")
          }
          ind = ind+1
          models[[ind]] = m
          modelnames = c(modelnames,paste0("SA",as.character(i),".",mname,".",dname))
        }
      }
    }
  }
}

names(models) = modelnames
save(modelnames, file = "modelnames.RData")

save(models,file ="unfitted_models.RData")


