setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##Step 1: download files from worldclim
library(raster)
library(sp)
library(dplyr)
library(Directional)
library(ggplot2)
library(raster)
library(viridis)
library(rasterVis)
library(httr)



##library(ShiftRotate)  PLEASE LOAD THE  FUNCTIONS FROM THE SHIFT & ROTATE SCRIPT DIRECTLY
#THE PACKAGE IS AVAILABLE ON GITHUB BUT THE REPO IS NOT ANONYMOUS

#if already dowloaded the file set to TRUE
set.seed(1)
download.files = T

#download the elevation file
if(download.files){
  url <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_elev.zip"
  destfile <- "elv.zip"
  download.file(url, destfile)
  outDir<-"elv"
  unzip(destfile,exdir=outDir)
}

#process elevation for original study zone
s=raster('elv/wc2.1_10m_elev.tif')
k <- crop(s, extent(8, 15, -6, 6))
o <- aggregate(k, fact = 3)
vartable <- as.data.frame(flip(o, direction = 'y'), xy = TRUE)
vartable$elv <- vartable[,3]
vartable[,3] <- NULL

#download other 19 bioclimatic variables from worldclim
if(download.files){
  url <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"
  destfile <- "bio.zip"
  download.file(url, destfile)
  outDir<-"bio"
  unzip(destfile,exdir=outDir)
}


#process the 5 chosen variables for Earth
varlist<- c(4, 5, 6, 7, 12)
linklist <- c('bio/wc2.1_10m_bio_4.tif',
              'bio/wc2.1_10m_bio_5.tif', 
              'bio/wc2.1_10m_bio_6.tif', 
              'bio/wc2.1_10m_bio_7.tif', 
              'bio/wc2.1_10m_bio_12.tif')


q=1
for( q in 1:length(linklist)){
  s=raster(linklist[[q]])
  k <- crop(s, extent(8, 15, -6, 6))
  o <- aggregate(k, fact = 3)
  v <- as.data.frame(flip(o, direction = 'y'), xy = TRUE)
  d <- data.frame(v[,3])
  colnames(d) <- paste('BIO',varlist[q], sep="")
  vartable <- cbind(vartable, d)
}



#make the raster layers for later 
ELV<- aggregate(raster('elv/wc2.1_10m_elev.tif'), fact = 3)
BIO4R<- aggregate(raster('bio/wc2.1_10m_bio_4.tif'), fact = 3)
BIO5R<- aggregate(raster('bio/wc2.1_10m_bio_5.tif'), fact = 3)
BIO6R<- aggregate(raster('bio/wc2.1_10m_bio_6.tif'), fact = 3)
BIO7R<- aggregate(raster('bio/wc2.1_10m_bio_7.tif'), fact = 3)
BIO12R<- aggregate(raster('bio/wc2.1_10m_bio_12.tif'), fact = 3)

#process elevation for Earth
s=raster('elv/wc2.1_10m_elev.tif')
o <- aggregate(s, fact = 3)
fulltable <- as.data.frame(flip(o, direction = 'y'), xy = TRUE)
fulltable$elv <- fulltable[,3]
fulltable[,3] <- NULL

#process the 5 chosen variables for Earth
varlist<- c(4, 5, 6, 7, 12)
linklist <- c('bio/wc2.1_10m_bio_4.tif',
              'bio/wc2.1_10m_bio_5.tif', 
              'bio/wc2.1_10m_bio_6.tif', 
              'bio/wc2.1_10m_bio_7.tif', 
              'bio/wc2.1_10m_bio_12.tif')

for( q in 1:length(linklist)){
  s=raster(linklist[[q]])
  o <- aggregate(s, fact = 3)
  v <- as.data.frame(flip(o, direction = 'y'), xy = TRUE)
  d <- data.frame(v[,3])
  colnames(d) <- paste('BIO',varlist[q], sep="")
  fulltable <- cbind(fulltable, d)
}


############ Shift and Rotate ###############

##quick visualizaiton check
coord = data.frame(vartable[, 2], vartable[, 1])
Lat.range = c(-90, 90)
Lon.range = c(-180, 180)
rotation = Lat.shift = Lon.shift = T
mirror = "r"

 plot(
   coord[, 2:1],
   asp = 1,
   xlim = c(-180, 180),
   ylim = c(-90, 90)
 )
 lines(cbind(
   c(Lon.range[1], Lon.range[2], Lon.range[2], Lon.range[1], Lon.range[1]),
   c(Lat.range[1], Lat.range[1], Lat.range[2], Lat.range[2], Lat.range[1])
 ))

 sphereplot(LonLat2XYZ(coord), col = 1) #projection on a sphere

### get the 100 datasets

out <- c()

for (i in 1:2000) {
  cord <-
    shiftrotGlobe(
      coord,
      Lat.range = Lat.range,
      Lon.range = Lon.range,
      Lat.shift = Lat.shift,
      Lon.shift = Lon.shift,
      rotation = rotation,
      mirror = mirror
    )
  
  #save extracts from polygons
  outE <- extract(ELV,   cord, df = TRUE, cellnumbers = TRUE)
  out4 <- extract(BIO4R, cord, df = TRUE, cellnumbers = TRUE)
  out5 <- extract(BIO5R, cord, df = TRUE, cellnumbers = TRUE)
  out6 <- extract(BIO6R, cord, df = TRUE, cellnumbers = TRUE)
  out7 <- extract(BIO7R, cord, df = TRUE, cellnumbers = TRUE)
  out12 <- extract(BIO12R, cord, df = TRUE, cellnumbers = TRUE)
  
  #save coordinates
  fx <- vartable$x
  ty <- vartable$y
  
  outp <-
    data.frame(
      x = fx,
      y = ty,
      elv = outE[, 3],
      BIO4 = out4[, 3],
      BIO5 = out5[, 3],
      BIO6 = out6[, 3],
      BIO7 = out7[, 3],
      BIO12 = out12[, 3]
    )
  
  #deleting data to match study zone with coastline
  missing <- is.na(rowSums(vartable))
  outp <- outp[!missing,]
  
  #saving output
  out[[i]] <- outp
}

# #trimming output down to 100 full datasets with 270 sample units
temp <- c()
for (i in 1:length(out)) {
  temp <-  append(temp, any(is.na(out[[i]])))
}


final <- out[!temp]
env_variables <- final[1:100]

#remove na values
varb <- !is.na(vartable[3])
vartable <- vartable[varb, ]

for (i in 1:length(env_variables)) {
  env_variables [[i]] <- (env_variables [[i]][varb, ])
}


#save the output in RData files
save(vartable, file = "environmental.data.original.RData")
save(env_variables , file = "environmental.data.shifed.and.rotated.RData")
