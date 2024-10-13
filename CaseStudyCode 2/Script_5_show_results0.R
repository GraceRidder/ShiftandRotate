setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##helpful function 
gelStats <- function(BetaGelman) {

  if (length(BetaGelman) > 400) {
    glmns <- c()
    x <- 7
    repeat {
      glmns <- append(glmns, mean(BetaGelman[, 1][x - 6:x]))
      x <- x + 7
      if (x > 700) {
        break
      }
    }
  } else {
    glmns <- c()
    x <- 2
    repeat {
      glmns <- append(glmns, mean(BetaGelman[, 1][x - 1:x]))
      x <- x + 2
      if (x > 200) {
        break
      }
    }
  }


  non_converge.sp <- c()
  converge <- c()
  for (i in 1:length(glmns)){
    if (glmns[i] > 1.1){
      non_converge.sp <- append(non_converge.sp, paste0("SP_",as.character(i)))
    } else {
      converge <- append(converge, paste0("SP_",as.character(i)))
    }
  }

  results=list(whichNotconverged = non_converge.sp,
               Notconverged = length(non_converge.sp),
               converge = converge)
  return(results)

}

#####  Define the model directory with the results of Script_4_fit_and_evaluate_models0.R
model.directory  <- "~/Desktop/HMSC_scripts/"
setwd(model.directory)

#organizing files
ts = 2
i = 1
samples = 250
thin = 1

fileindex <- c()
for (i in 1:4800) {
  filename = file.path(
    model.directory,
    paste0(
      "results_",
      as.character(ts),
      "/",
      "model_",
      as.character(i),
      "_samples_",
      as.character(samples),
      "_thin_",
      as.character(thin),
      ".RData"
    )
  )
  fileindex <- append(fileindex, filename)
}

load("modelnames.RData")

x <- c()
a <- 1:4800
for (i in 1:4) {
  b <- a[seq(i, length(a), 12)]
  x <- c(x, b)
}

X <- sort(x)
OG16 <- fileindex[X]
OGnames <- modelnames[X]

x <- c()
for (i in 5:8) {
  b <- a[seq(i, length(a), 12)]
  x <- c(x, b)
}
X <- sort(x)
Within <- fileindex[X]
W.names <- modelnames[X]

x <- c()
for (i in 9:12) {
  b <- a[seq(i, length(a), 12)]
  x <- c(x, b)
}

X <- sort(x)
Across <- fileindex[X]
A.names <- modelnames[X]

###CHOOSE ONE FOR AUC GRAPHS
models <- Across
names <- A.names

models <- Within
names <- W.names

models <- OG16
names <- OGnames

AUClist <- c()
Gellist <-  list()
OGBmat <- list()
OGM <- list()
OGsup <- list()

for (i in 1:1600) {
   load(models[i])
    AUClist <- rbind(AUClist, AUC.tot)
    Gellist[[i]] <- BetaGelman
    OGM[[i]] <- postBeta$mean
    OGsup[[i]] <-  postBeta$support

}
colnames(AUClist) <- c("MF", "CV", "CVS", "CVE")

####visualizing AUC values
library(ggplot2)
x <- c("M.E6", "M.E6.S", "M.E1", "M.E1.S")
y <- c("D.E6.P", "D.E6.S.P", "D.E1.P", "D.E1.S.P")
data <- expand.grid(Modeltype=x, Datatype=y)

st <-names[1:16]
new_str <- gsub("SA1.","",st, fixed = T)
list <- new_str

evens <- function(x) subset(x, x %% 2 == 0)
e <- evens(1:200)
ee <- e-1

AUC <- c()
for (i in 1:16) {
  c <- grep(list[i], names, fixed = T)
  if (length(c)>100){
    c <- c[ee]}
  m <- AUClist[, 1][c]  ######  change number from 1-4 to generate graphs with different AUC value types. 
  AUC <- append(AUC, mean(m)) #### order is 1:AUC, 2:cross validation 3:cross validation with space, 4:cross validation with environment 
}







ggplot(data, aes(Datatype, Modeltype, fill= AUC)) +
  geom_tile() +
  geom_text(aes(label = round(AUC, digits = 3))) +
  scale_fill_gradient(low = "yellow", high = "red", limits = c(0, 1)) +
  ggtitle("f)") +  ###change name to fit model type 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()  )



#true parameters
#load("datasets.RData")

trove <- unlist(true.parameters, recursive = F)

tops <- c()
for (k in 1:length(trove)) {
  for (i in 1:4) {
    tops <- append(tops, trove[k])
  }
}

topper <- c()
for (k in 1:length(trove)) {
  for (i in 1:4) {
    topper <- append(topper, trove[k])
  }
}

xv <- c()
x <- -1
v <- 0

while(v < 1600){
  x <- x + 4
  v <- v + 4
  xv <- append(xv, x)
  xv <- append(xv, v)
}

for (o in xv){
  topper[[o]] <- topper[[o]][-c(3:7),]
}

r <- -7
e <- -6
g <- -3
w <- -2

regw <- c()
while(w < 1598){
  r <- r + 16
  e <- e + 16
  g <- g + 16
  w <- w + 16
  regw <- append(regw, r)
  regw <- append(regw, e)
  regw <- append(regw, g)
  regw <- append(regw, w)
}

top <- topper

for (o in regw){{
  topper[[o]][3:7,]  <- 0
  }
}


#turning support binary
temp <- OGsup
for (o in 1:length(temp)){
  for (k in 1:length(temp[[o]])){
    if (OGsup[[o]][k] > .95 | OGsup[[o]][k] < .05 ){
      temp[[o]][k] <- 1
    } else {
      temp[[o]][k] <- 0
    }
  }
}

for (o in 1:length(topper)){
  temp[[o]] <- temp[[o]][-c(1),]
  topper[[o]] <- topper[[o]][-c(1), ]
  OGM[[o]] <-OGM[[o]][-c(1),]
  
}

###estimating correct, incorrect and not statistically supported
correctTlist <- c()
incorrectTlist <- c()
unsupportedlist <- c()
correctZlist <- c()
incorrectZlist <- c()
nosupZlist <- c()
truezeroslist<- c()
supportedlist <- c()

correctXlist <- c()
incorrectXlist <- c()
nosupXlist <- c()
nonzeroslist <- c()

###funneling down to specific model-data combination
for (j in 1:16) {
  
  c <- grep(list[j], names)
  if (length(c)>100){
    c <- c[ee]
  }
  
  support <- temp[c]
  estimated <- OGM[c]
  trueB <- topper[c]

  ### etc
  for (n in 1:length(trueB)){
    correctT = 0
    incorrectT =0
    unsupported=0
    nosupZ=0
    incorrectZ=0
    correctZ=0
    nosupX=0
    incorrectX=0
    correctX=0
    nonzeros = 0
    truezeros = 0
    supported = 0

    for (i in 1:length(trueB[[n]])) {

      ### All values
      #correct
      if (sign(trueB[[n]][i]) == sign(estimated[[n]][i]) && support[[n]][i] == 1)  {
        correctT <- correctT + 1
      }
      #incorrect
      if (sign(trueB[[n]][i]) != sign(estimated[[n]][i]) && support[[n]][i] == 1) {
        incorrectT <- incorrectT + 1
      }
      #not supported
      if (support[[n]][i] == 0) {
        unsupported <- unsupported + 1
      }

      ### True value zero
      #No support
      if (support[[n]][i] == 0 && sign(trueB[[n]][i]) == 0) {
        nosupZ <- nosupZ + 1
      }
      #support and estimate does not equal zero
      if (support[[n]][i] == 1 && sign(trueB[[n]][i]) == 0 && sign(estimated[[n]][i]) != 0) {
        incorrectZ <- incorrectZ + 1
      }
      #support and estimate does equal zero
      if (support[[n]][i] == 1 && sign(trueB[[n]][i]) == 0 && sign(estimated[[n]][i]) == 0) {
        correctZ <- correctZ + 1
      }

      ### True value is nonzero ###
      #No support
      if (support[[n]][i] == 0 && sign(trueB[[n]][i]) != 0) {
        nosupX <- nosupX + 1
      }
      #support and estimate does not equal zero
      if (support[[n]][i] == 1 && sign(trueB[[n]][i]) != 0 && sign(estimated[[n]][i]) != sign(trueB[[n]][i])) {
        incorrectX <- incorrectX + 1
      }
      #support and estimate does equal zero
      if (support[[n]][i] == 1 && sign(trueB[[n]][i]) != 0 && sign(estimated[[n]][i]) == sign(trueB[[n]][i])) {
        correctX <- correctX + 1
      }

      if (support[[n]][i] == 1){
        supported <- supported + 1
      }
      if (sign(trueB[[n]][i]) == 0){
        truezeros <- truezeros + 1
      }
      if (sign(trueB[[n]][i]) != 0){
        nonzeros <- nonzeros + 1
      }
    }
    
   ###these are total 
    correctTlist <- append(correctTlist, correctT)
    incorrectTlist <- append(incorrectTlist, incorrectT)
    unsupportedlist <- append(unsupportedlist, unsupported)

    ##these are the zeros (compare to true zeros)
    correctZlist <- append(correctZlist, correctZ)
    incorrectZlist <- append(incorrectZlist, incorrectZ)
    nosupZlist <- append(nosupZlist, nosupZ)

    #these are non zeros, compare to nonzeros 
    correctXlist <- append(correctXlist, correctX)
    incorrectXlist <- append(incorrectXlist, incorrectX)
    nosupXlist <- append(nosupXlist, nosupX)
    
    
    nonzeroslist <- append(nonzeroslist, nonzeros)
    truezeroslist <- append(truezeroslist, truezeros)
    supportedlist <- append(supportedlist, supported)
  }
    }

## CHANGE THIS BETWEEN LISTS FOR DIFFERENT GRAPHS
true <- (correctXlist)
inco <- (incorrectXlist)
nosu <- (nosupXlist)

tot <- true + inco + nosu

## graphs
dd <- colMeans(matrix(true,100))  #### change from true, inco, or nosu.
zz <- colMeans(matrix(tot,100)) 
percent <- dd/zz


ggplot(data, aes(Datatype, Modeltype, fill= percent)) +
  geom_tile() +
  geom_text(aes(label = round(percent, digits = 3))) +
  scale_fill_gradient(low = "yellow", high = "red", limits = c(0, 1)) +
  ggtitle("correct zeros") ### CHANGE NAME 

