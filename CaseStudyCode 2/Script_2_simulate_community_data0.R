library(MASS)
set.seed(1)
load("environmental.data.original.RData")
load("environmental.data.shifed.and.rotated.RData") # named env_variables
env_variables[[1]] = vartable #set the African landscape as the first study area
nd = length(env_variables)

data.sets = list()
true.parameters = list()
environmental.covariates = list()
spatial.coordinates = list()
for(i in 1:nd){
  a = env_variables[[i]]
  xy = data.frame(x=a$x,y=a$y)
  di = as.matrix(dist(xy))
  X = a
  X$x = NULL
  X$y = NULL
  X = scale(X)
  n = dim(X)[1]
  betas = list()
  ind = 0
  data.set = list()
  for(env in c(".E6",".E1")){
    for(spa in c("",".S")){
      ind = ind + 1
      if(env=="E1"){
        beta = matrix(rnorm(100*1,sd = 1),nrow = 1)
        X1 = X[,1]
      } else {
        beta = matrix(rnorm(100*6,sd = sqrt(1/6)),nrow = 6)
        X1 = X
      }
      X1 = cbind(1,X1)
      beta = rbind(matrix(rnorm(100*1,sd = 1),nrow = 1),beta)
      LF = X1%*%beta
      if(spa==".S"){
        LR = mvrnorm(n = 1, mu=rep(0,n), Sigma=exp(-di/3))
        L = LF + LR
      } else {
        L = LF
      }
      eps = matrix(rnorm(100*n),nrow = n)
      Y = 1*((L+eps)>0)
      data.set[[ind]] = Y
      names(data.set)[ind] = paste0("D",env,spa)
      betas[[ind]] = beta
      names(betas)[ind] = paste0("D",env,spa)
    }
  }
  data.sets[[i]] = data.set
  true.parameters[[i]] = betas
  spatial.coordinates[[i]] = xy
  environmental.covariates[[i]] = X
}
save(data.sets,true.parameters,spatial.coordinates,environmental.covariates,file = "datasets.RData")
