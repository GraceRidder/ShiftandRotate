###################################################################################################################
# Code illustrating the usage and possible applications of the Shift & Rotate functions

# Written by Olivier J. Hardy for the article by Ridder, G., Hardy, O. J., Ovaskainen, O., entitled 
# "Generating spatially Realistic Environmental Null Models with the Shift-&-Rotate Approach Helps 
# Evaluate False Positives in Species Distribution Modeling", Methods in Ecology and Evolution (submitted)
##########################################################################################################################


source("Shift&RotFunctions.R")

#############################################################################################
# Test functions: 1 shiftrotXY(), 2 shiftrotGlobe(), 3 rrotGlobe()
# to visualize how they move spatial coordinates

library(Directional) #to visualize points on a sphere using sphereplot()
library(geodata)
library(sdm)


# define coordinates of 7 points forming a "half-arrow"
coord=data.frame(lon=c(-80,-80,-80,-80,-80,-80,-75),lat=c(1,6,11,16,21,26,6))
# define the window within which coord must be moved
Lat.range = c(-20,80)
Lon.range = c(-90,-10)

# set functions parameters for random rotation and longitudinal/latitudinal shifts, and randomly using a mirror image
rotation = Lat.shift = Lon.shift = T
mirror = "r"

# plot original coord and window
plot(coord, asp=1, xlim=c(-180,180), ylim=c(-90,90))
lines( cbind(c(Lon.range[1],Lon.range[2],Lon.range[2],Lon.range[1],Lon.range[1]), c(Lat.range[1],Lat.range[1],Lat.range[2],Lat.range[2],Lat.range[1])))
sphereplot(LonLat2XYZ(coord), col=1) #projection on a sphere that should appear as a 'RGL' window where the sphere can be rotated using the mouse 

# plot 6 times randomized coord using one of the three possible functions:
# 1 shiftrotXY(), 2 shiftrotGlobe(), 3 rrotGlobe()
fun = 2 #try 1, 2 or 3
for(i in 1:6){ 
  if(fun==1) newcoord = shiftrotXY(coord[,1:2], X.range=Lon.range, Y.range=Lat.range, X.shift=Lon.shift, Y.shift=Lat.shift, rotation=rotation, mirror=mirror)[,1:2]
  if(fun==2) newcoord = shiftrotGlobe(coord, Lon.range=Lon.range, Lat.range=Lat.range, Lon.shift=Lon.shift, Lat.shift=Lat.shift, rotation=rotation, mirror=mirror)
  if(fun==3) newcoord = rrotGlobe(coord, Lon.range=Lon.range, Lat.range=Lat.range, mirror=mirror)
  points(newcoord, col=1+i)
  sphereplot(LonLat2XYZ(newcoord), col=1+i)
}

# one can try also with
Lat.range = c(-90,90); Lon.range = c(-180,180) #unlimited lat/lon ranges
Lat.range = c(30,70); Lon.range = c(-90,10)    #window outside original lat/lon ranges
Lat.range = c(-70,0); Lon.range = c(160,-160)  #longitude range crossing the 180°/-180° limit
Lat.range = c(0,15); Lon.range = c(0,15)  #too small window



#############################################################################################
# Example of usage of shiftrotGlobe() to rotate or shift longitudinally and latitudinally
# a set of points by systematic angles

coord=data.frame(long=c(-80,-80,-80,-80,-80,-80,-75), lat=c(1,6,11,16,21,26,6))
Lon.range = c(-180,180)
Lat.range = c(-90,90)
plot(coord, asp=1, xlim=Lon.range, ylim=Lat.range)
lines( cbind(c(Lon.range[1],Lon.range[2],Lon.range[2],Lon.range[1],Lon.range[1]), c(Lat.range[1],Lat.range[1],Lat.range[2],Lat.range[2],Lat.range[1])))

# set initial function parameters to have no rotation and no longitudinal/latitudinal shifts
rotation = Lat.shift = Lon.shift = mirror = F
# Plot coord after systematic 45° rotations
for(i in 1:7){ 
  rotation = i*45
  newcoord = shiftrotGlobe(coord, Lon.range=Lon.range, Lat.range=Lat.range, Lon.shift=Lon.shift, Lat.shift=Lat.shift, rotation=rotation, mirror=mirror)
  points(newcoord, col=1+i)
}
# Plot coord after systematic 20° shifts in longitude 
coord2 = newcoord
rotation = F
for(i in 1:17){ 
  Lon.shift = i*20
  newcoord = shiftrotGlobe(coord2, Lon.range=Lon.range, Lat.range=Lat.range, Lon.shift=Lon.shift, Lat.shift=Lat.shift, rotation=rotation, mirror=mirror)
  points(newcoord, col=1+i)
}
# Plot coord after systematic 20° shifts in latitude 
Lon.shift = F
for(i in 1:17){ 
  Lat.shift = i*20 
  newcoord = shiftrotGlobe(coord2, Lon.range=Lon.range, Lat.range=Lat.range, Lon.shift=Lon.shift, Lat.shift=Lat.shift, rotation=rotation, mirror=mirror)
  points(newcoord, col=1+i)
}






#############################################################################################
# Example of usage of shiftrotGlobe() to generate null datasets for SDM application.
# The script uses GBIF to extract occurrence data of a species, WORLDCLIM to extract climatic 
# layers at low resolution (10min degrees) as predictors at the world scale, and the R package 
# 'sdm' to model the species distribution using different algorithms (glm, gam, random forest). 
# Statistics summarizing model performance (AUC and COR in the present case) computed 
# for the real dataset and 100 null datasets obtained by shiftrotGlobe() allow assessing .
# the risk over model over-fitting leading to an excess of false positives. 
# The SDM procedure implemented hereafter is basic and does not pretend to be optimal for 
# modeling species distribution, the goal being solely to illustrate how to apply the 
# Shift&Rotate approach.  

library(terra)
library(sp)
library(geosphere)

#download occurrence data from GBIF
#coordSpGBIF = sp_occurrence("Pericopsis", "elata")
coordSpGBIF = sp_occurrence("Baillonella", "toxisperma")

plot(coordSpGBIF[,c("lon","lat")], asp=1)
#remove points outside realistic natural range (to adjust according to the species) and keep only lon-lat coord
coordSp = coordSpGBIF[coordSpGBIF$lat>-5 & coordSpGBIF$lat<10 & coordSpGBIF$lon>8 & coordSpGBIF$lon<20 & !is.na(coordSpGBIF$lon), c("lon","lat")]
plot(coordSp, asp=1)

#generate 500 random background points within a window 5° broader than the lon-lat ranges of the species
coordBg = data.frame( lon = runif(500, min(coordSp$lon)-5, max(coordSp$lon)+5 ), 
                      lat = runif(500, min(coordSp$lat)-5, max(coordSp$lat)+5 )   )
#remove background points <30 km from an occurrence point
distBgSp = distm(coordBg, coordSp)
mindistBgSp = apply(distBgSp, 1, FUN=min)
coordBg = coordBg[mindistBgSp > 50000, ]

#merge coordSp and coordBg and add a column for presence/absence
coord = rbind(coordSp, coordBg)
coord$sp = c(rep(1, nrow(coordSp)), rep(0, nrow(coordBg)))

#download worldclim data at lowest resolution
WC2 = worldclim_global(var="bio", res=10, path = getwd() )
# load WorldClim layers
lst <- list.files(path="./climate/wc2.1_10m", pattern='tif', full.names = T)
lst = lst[c(1,12:19,2:11)] #Bioclim data
pred <- rast(lst[c(1,2,3,12,14,15)]) #BioClim n°
#pred <- rast(lst) #all Biolayers
plot(pred)

#crop predictors to coord and show occurrences and background points 
predc= crop(pred, ext(min(coord$lon), max(coord$lon), min(coord$lat), max(coord$lat)))
plot(predc)
plot(predc[[1]]) #plot annual precipitation
points(coord[,1:2], pch= c(3,1)[coord$sp+1], col=c("red", "blue")[coord$sp+1] )

#remove coordinates undefined for predictors
pred_coord = extract(pred, coord[,1:2])
Nmiss = rowSums(is.na(pred_coord))
coord = coord[Nmiss == 0,]

#create sdm object that will be used to run the SDM models
sdmD <- sdmData(~.+coord(lon+lat), train=coord, predictors=pred)

# fit models with three methods (glm = generalized linear model, gam = generalized additive model, rf = random forest)
m <- sdm(data=sdmD, methods=c('glm','gam', 'rf'))

m   # summary statistics of the different models using the real predictors

# Generate a null dataset with the Shift & Rotate approach
coordR = coord #coordR will store the set of randomized coordinates while keeping the relative positions of data points
Nmiss=1
while( Nmiss ){ # while condition to repeat the shift&rotate moves until all points fall where climate predictors are defined 
  coordR[,1:2] = shiftrotGlobe(coord[,1:2]) #by default the set of points are randomly moved around the globe
  pred_R = extract(pred, coordR[,1:2])  #extract climate variables (predictors) for the set of points in coordR
  Nmiss = sum(is.na(pred_R)) #number of missing data for predictors after moving the points 
}
sdmDR <- sdmData(~.+coord(lon+lat), train=coordR, predictors=pred)
plot(sdmDR, asp=1)

mR <- sdm(data=sdmDR, methods=c('glm','gam', 'rf'))
mR    # summary statistics of the different models using 'fake' but realistic predictors

# Extract AUC and COR stats for the model with real data
AUCobs = data.frame( glm = m@models$sp$glm$'1'@evaluation$training@statistics$AUC,
                     gam = m@models$sp$gam$'2'@evaluation$training@statistics$AUC,
                     rf = m@models$sp$rf$'3'@evaluation$training@statistics$AUC       )
CORobs = data.frame( glm = m@models$sp$glm$'1'@evaluation$training@statistics$COR[1],
                     gam = m@models$sp$gam$'2'@evaluation$training@statistics$COR[1],
                     rf = m@models$sp$rf$'3'@evaluation$training@statistics$COR[1]    )

# Make 100 S&R null models to assess the null distributions of AUC and COR
AUCr = CORr = as.data.frame(matrix(nrow=100, ncol=3))
colnames(AUCr) = colnames(CORr) = c("glm", "gam", "rf")
# Define the window within which the Shift&Rotate ()
Lat.range = c(-90,90); Lon.range = c(-180,180) #unlimited lat/lon ranges
Lat.range = c(-23,23); Lon.range = c(-10,40) #limited to tropical Africa

for(r in 1:100){
  cat(r, sep=" ") #counter to check progress in analyzing null datasets
  coordR = coord
  Nmiss=1
  while( Nmiss ){
    coordR[,1:2] = shiftrotGlobe(coord[,1:2], Lon.range=Lon.range, Lat.range=Lat.range)
    pred_R = extract(pred, coordR[,1:2])
    Nmiss = sum(is.na(pred_R))
  }
  sdmDR <- sdmData(~.+coord(lon+lat), train=coordR, predictors=pred)
  mR <- sdm(data=sdmDR, methods=c('glm','gam','rf'))
  # copy statistics after verifying if the model worked
  if(length(mR@models$sp$glm$'1'@evaluation)){
    AUCr[r,1] = mR@models$sp$glm$'1'@evaluation$training@statistics$AUC
    CORr[r,1] = mR@models$sp$glm$'1'@evaluation$training@statistics$COR[1]
  }
  if(length(mR@models$sp$gam$'2'@evaluation)){ 
    AUCr[r,2] = mR@models$sp$gam$'2'@evaluation$training@statistics$AUC
    CORr[r,2] = mR@models$sp$gam$'2'@evaluation$training@statistics$COR[1]
  }
  if(length(mR@models$sp$rf$'3'@evaluation)){
    AUCr[r,3] = mR@models$sp$rf$'3'@evaluation$training@statistics$AUC
    CORr[r,3] = mR@models$sp$rf$'3'@evaluation$training@statistics$COR[1]
  }
}

# Compute proportion of null models with AUC or COR at least as high as with the observed data
Pval = sum(AUCr$glm >= AUCobs$glm, na.rm=T) / sum(!is.na(AUCr$glm))
hist(AUCr$glm, main= paste("AUCr - glm  Pval=", Pval) )
abline(v = AUCobs$glm)

Pval = sum(AUCr$gam >= AUCobs$gam, na.rm=T) / sum(!is.na(AUCr$gam))
hist(AUCr$gam, main= paste("AUCr - gam  Pval=", Pval) )
abline(v = AUCobs$gam)
     
Pval = sum(AUCr$rf >= AUCobs$rf, na.rm=T) / sum(!is.na(AUCr$rf))
hist(AUCr$rf, main= paste("AUCr - rf  Pval=", Pval) )
abline(v = AUCobs$rf)

Pval = sum(CORr$glm >= CORobs$glm, na.rm=T) / sum(!is.na(CORr$glm))
hist(CORr$glm, main= paste("CORr - glm  Pval=", Pval) )
abline(v = CORobs$glm)

Pval = sum(CORr$gam >= CORobs$gam, na.rm=T) / sum(!is.na(CORr$gam))
hist(CORr$gam, main= paste("CORr - gam  Pval=", Pval) )
abline(v = CORobs$gam)

Pval = sum(CORr$rf >= CORobs$rf, na.rm=T) / sum(!is.na(CORr$rf))
hist(CORr$rf, main= paste("CORr - rf  Pval=", Pval) )
abline(v = CORobs$rf)


# Tentative interpretation:
# gam and rf usually lead to higher performance (higher AUC and COR) than glm models.
# However, the AUC of gam and rf can be very high even using fake but realistic predictors, 
# with sometimes AUC = 1 for all S&R null datasets, as for the rf model.
# When comparing the AUC or COR observed with the real data with their corresponding 
# distributions for the S&R null datasets, the Pvalue (i.e. proportion of AUC or COR obtained
# under the S&R null datasets >=  AUC or COR observed with the real dataset) # can be lower 
# for glm than gam or rf, indicating a higher power of glm to identify causal predictors.
# Hence, the apparent higher performance of gam or rf models based on their high AUc or COR 
# result essentially from over-fitting, while they can potentially be less good than 
# simple glm to correctly infer causal predictors.  
