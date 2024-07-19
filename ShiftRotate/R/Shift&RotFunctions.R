#########################################################################################################################
# Functions to apply (random) rotation and/or shift of a set of spatial coordinates preserving their relative
# spatial distribution (i.e. distance matrix) on a rectangular area or on a sphere. The resulting spatial
# coordinates can be constrained to fit within a defined range of X-Y or longitude-latitude coordinates.
# Optionally, one coordinate can also be reversed to obtain the mirror image.
#
# Three main functions ('shiftrotXY', 'shiftrotGlobe', 'rrotGlobe') aim to help programming null models of realistic 
# spatial predictors (e.g. climate data) defined over a continuous area (raster of environmental layers, globe). 
# They return a set of moved spatial coordinates.
#
# Three additional functions ('LonLat2XYZ', 'XYZ2LonLat', 'RotLonLatAxis') are used by the main functions 
# to convert spatial coordinates on a sphere between Cartesian (X, Y, Z) and polar (longitude, latitude) 
# coordinate systems, and to rotate a set of coordinates on a sphere around a single axis.
# 
# Written by Olivier J. Hardy for the article by Ridder, G., Hardy, O. J., Ovaskainen, O., entitled 
# "Generating spatially Realistic Environmental Null Models with the Shift-&-Rotate Approach Helps 
# Evaluate False Positives in Species Distribution Modeling", Methods in Ecology and Evolution (submitted)
##########################################################################################################################


##########################################################################################################################
# LonLat2XYZ(coord) converts longitude, latitude, coordinates in degrees into Cartesian XYZ coordinates 
# (where X points to 0° latitude and 0° longitude, Y points to 0° latitude and +90° longitude, Z points to 90° latitude). 
LonLat2XYZ = function(coord){ 
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  coord = coord/180
  a1 = sinpi(1/2 - coord[, 2])
  XYZ = cbind(a1 * cospi(coord[, 1]), a1 * sinpi(coord[, 1]), cospi(1/2-coord[, 2]))
  colnames(XYZ) = c("X", "Y", "Z")
  XYZ
}

##########################################################################################################################
# XYZ2LonLat(XYZ) converts Cartesian XYZ coordinates on a sphere (satisfying X^2 + Y^2 + Z^2 = 1) into 
# longitude (first column, ]-180°, 180°]) and latitude (second column, [-90°, 90°]) coordinates in degrees. 
XYZ2LonLat = function(XYZ){ 
  XYZ = as.matrix(XYZ)
  if (ncol(XYZ) == 1) XYZ = t(XYZ)
  Lat = 90 - (acos(XYZ[,3]) * 180/pi)
  Lon = ( (atan(XYZ[, 2]/XYZ[, 1]) + pi*I(XYZ[,1] < 0))%%(2*pi)) * 180/pi 
  Lon[XYZ[,3]==1] = 0
  Lon[Lon>180] = Lon[Lon>180] - 360
  cbind(Lon,Lat)
}

##########################################################################################################################
# RotLonLatAxis(coord, axis, angle) rotates a set of longitude, latitude coordinates around a given axis by a given angle.  
RotLonLatAxis = function(coord, axis, angle){
  R_n_a = function(a, n1, n2, n3)  rbind( c( cos(a)+n1^2*(1-cos(a)), n1*n2*(1-cos(a))-n3*sin(a), n1*n3*(1-cos(a))+n2*sin(a)),
                                          c( n1*n2*(1-cos(a))+n3*sin(a), cos(a)+n2^2*(1-cos(a)), n2*n3*(1-cos(a))-n1*sin(a)),
                                          c( n1*n3*(1-cos(a))-n2*sin(a), n2*n3*(1-cos(a))+n1*sin(a), cos(a)+n3^2*(1-cos(a)))  ) 
  XYZ = LonLat2XYZ(coord)
  xyz = LonLat2XYZ(axis)
  Rna = R_n_a(-1*angle/180*pi, xyz[1], xyz[2], xyz[3])
  RXYZ = XYZ %*% Rna
  Rcoord = XYZ2LonLat(RXYZ)
  Rcoord
}


##########################################################################################################################
# shiftrotXY() performs (random) rotation and/or shift on a rectangular area where position is defined by 
# Cartesian X, Y coordinates. It can be applied using longitude and latitude but will not correct for Earth curvature, 
# so that it is not recommended outside the inter-tropical zone and for an area spanning over several tens of degrees 
# longitude or latitude.  

shiftrotXY = function(coord, X.range=c(-180,180), Y.range=c(-30,30), rotation=T, X.shift = T, Y.shift = T, 
                      mirror = c("random", "no", "yes"), verbose = T, MAXtrial = 1000 ) 
{
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  Xp = coord[,1]
  Yp = coord[,2]
   #if no range is given, assume coordinates are degrees longitude (X) and latitude (Y). The algorithm is not optimal near the poles due to extreme area deformation

  if(verbose){
    #check if points are within X-Y ranges (this does not impede the rest of the algorithm). 
    if((out=sum(Xp<X.range[1] | Xp>X.range[2] | Yp<Y.range[1] | Yp>Y.range[2])))
      print(paste("Warning:",out,"original coordinates fall outside the X and Y ranges of the window"))
    #check if window is large enough for free rotation, translation
    maxdistp = max(dist(coord[,1:2]))
    if( (X.range[2] - X.range[1]) < maxdistp | (Y.range[2] - Y.range[1]) < maxdistp )
      print("Warning: window size defined by X.range and Y.range is too small to allow free rotation. Consider using a larger window.")
  }
  #moved coordinates
  Xm = Xp
  Ym = Yp
  
  #Step 1. if requested, take mirror image and recenter on original coordinates
  mir=F
  if(mirror==T | substr(mirror[1],1,1)=="y" | substr(mirror[1],1,1)=="Y" | ((substr(mirror[1],1,1)=="r" | substr(mirror[1],1,1)=="R") & runif(1)<0.5)) mir=T 
  if(mir) Xm = -Xp + max(Xp) + min(Xp)
  
  #Step 2. apply rotation around centroid, checking that points can enter the window size
  if(rotation){
      if(!is.logical(rotation) & is.numeric(rotation)) angle1 = rotation
      trial=0
      repeat{
        trial = trial + 1
        if(is.logical(rotation)) angle1 = runif(1)*360
        X = Xm - mean(Xm); Y = Ym - mean(Ym)
        Xr = X*cos(angle1/180*pi) - Y*sin(angle1/180*pi)
        Yr = X*sin(angle1/180*pi) + Y*cos(angle1/180*pi)
        Xm = Xr + mean(Xm); Ym = Yr + mean(Ym)
        if( (max(Xm) - min(Xm)) <= (X.range[2] - X.range[1]) &  (max(Ym) - min(Ym)) <= (Y.range[2] - Y.range[1]) ) break
        if(trial == MAXtrial | is.numeric(rotation)) break
      }
  } 
    
  #Step 3. apply X and Y translation
  if(X.shift){
    if(!is.logical(X.shift) & is.numeric(X.shift)){
      Xm = Xm + X.shift
    }else if((max(Xm) - min(Xm)) <= (X.range[2] - X.range[1])){
      dX = runif(1, min = X.range[1]-min(Xm), max = X.range[2]-max(Xm) )
      Xm = Xm + dX
    }
  }
  if(Y.shift){
    if(!is.logical(Y.shift) & is.numeric(Y.shift)){
      Ym = Ym + Y.shift
    }else if((max(Ym) - min(Ym)) <= (Y.range[2] - Y.range[1])){
      dY = runif(1, min = Y.range[1]-min(Ym), max = Y.range[2]-max(Ym) )
      Ym = Ym + dY
    }
  }
  
  #report results
  Rcoord = cbind(Xm, Ym)
  if(trial == MAXtrial) Rcoord = NULL
  Rcoord
}


##########################################################################################################################
# shiftrotGlobe() performs (random) rotation and/or longitude and latitude shifts of spatial coordinates on a sphere, 
# possibly within a limited window. Adapted to situations where the set of coordinates to move spans over less than 90° 
# longitude and latitude. 

shiftrotGlobe = function(coord, Lon.range=c(-180,180), Lat.range=c(-90,90), rotation=T, Lon.shift=T, Lat.shift=T, 
                         mirror=c("random"), verbose=T, MAXtrial=1000) 
{
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  Lon = coord[,1]
  Lat = coord[,2] 
  Lon.range.orig = Lon.range
  #operate longitudinal shift so that Lon.range starts at 0 
  Lon = (Lon - Lon.range.orig[1] + 360)%%360
  Lon.range = (Lon.range - Lon.range[1] + 360)%%360
  if(Lon.range[2] == 0) Lon.range[2] = 360
  #check if points are within lon-lat ranges (this does not impede the rest of the algorithm). 
  if(verbose){
    if((out=sum(Lon<Lon.range[1] | Lon>Lon.range[2] | Lat<Lat.range[1] | Lat>Lat.range[2])))
      print(paste("Warning:",out,"original coordinates fall outside the latitudinal or longitudinal ranges"))
  }
  #initialize moved coordinates
  RLon = Lon
  RLat = Lat
  
  #Step 1: if requested, take mirror image and recenter on original coordinates
  mir=F
  if(mirror==T | substr(mirror,1,1)=="y" | substr(mirror,1,1)=="Y" | ((substr(mirror,1,1)=="r" | substr(mirror,1,1)=="R") & runif(1)<0.5)) mir=T 
  if(mir) RLon = -Lon + max(Lon) + min(Lon)
  
  Rcoord = cbind(RLon, RLat)
  
  trial = 0
  repeat{
    trial = trial+1

    #Step 2: apply rotation around axis pointing to centroid of data points, checking that points can enter the window size
    #axis1 = centroid axis (mean vector of coordinates, scaled to unit length)
    axis1 = XYZ2LonLat( colMeans(LonLat2XYZ(Rcoord)) / sqrt(sum(colMeans(LonLat2XYZ(Rcoord))^2)) )
    if(rotation){
      if(!is.logical(rotation) & is.numeric(rotation)){
        R2coord = RotLonLatAxis(Rcoord, axis1, rotation)
      }else{
        angle1 = runif(1)*360
        R2coord = RotLonLatAxis(Rcoord, axis1, angle1)
      } 
    } else R2coord = Rcoord
    
    #Step 3: apply latitudinal rotation along an equatorial vector perpendicular to the centroid vector
    if(Lat.shift){
      axis2 = c( (axis1[,1] - 90 + 360)%%360 , 0 ) #equatorial axis (lat==0) pointing perpendicular to centroid axis1 
      if(!is.logical(Lat.shift) & is.numeric(Lat.shift)){
        R3coord = RotLonLatAxis(R2coord, axis2, Lat.shift)
    #    R3coord[,2] = (R3coord[,2] + 90)%%180 - 90 #to limit latitudes within the [-90, 90] range
      }else{
        angle2 = -acos(runif(1, min=sin(Lat.range[1]/180*pi), max=sin(Lat.range[2]/180*pi)))*180/pi + 90 #random latitude for a point uniformly distributed on a sphere between 2 latitudes
        angle2 = angle2 - axis1[2] #rescale angle so that a change in latitude from the centroid latitude leads to a probability of final latitude ~ cos(latitude) 
        R3coord = RotLonLatAxis(R2coord, axis2, angle2)
      }
    } else R3coord = R2coord
   
    #Step 4: apply longitudinal translation, adding a random longitudinal shift
    R4coord = R3coord
    if(Lon.shift){
      if(!is.logical(Lon.shift) & is.numeric(Lon.shift)){
        R4coord[,1] = (R4coord[,1] + Lon.shift + 360)%%360
      }else{
        angle3 = 0
        if((Lon.range[2]-Lon.range[1])==360) angle3 = runif(1, min = 0, max=360)
        if((Lon.range[2]-Lon.range[1])<360 & (max(R3coord[,1])-min(R3coord[,1])) <= (Lon.range[2]-Lon.range[1]) ) 
          angle3 = runif(1, min = Lon.range[1]-min(R3coord[,1]), max = Lon.range[2]-max(R3coord[,1]))  
        R4coord[,1] = (R4coord[,1] + angle3 + 360)%%360
      } 
    }
    
    # Finally, check that moved coordinates fall in the window, otherwise try again, unless there is no random moves or MAXtrial is reached
    if(min(R4coord[,2])>=Lat.range[1] & max(R4coord[,2])<=Lat.range[2] &  
       min(R4coord[,1])>=Lon.range[1] & max(R4coord[,1])<=Lon.range[2]) break
    if(rotation!=T & Lat.shift!=T & Lon.shift!=T) break 
    if(trial == MAXtrial) break
  }
  
  if(verbose){
    if( min(R4coord[,2]) < Lat.range[1] | max(R4coord[,2]) > Lat.range[2] )
      print(paste("Warning: failed to move original coordinates into the latitudinal range after",trial,"trials."))
    if(  (Lon.range[2]-Lon.range[1])<360 & (min(R4coord[,1]) < Lon.range[1] | max(R4coord[,1]) > Lon.range[2]) )
      print(paste("Warning: failed to move original coordinates into the longitudinal range after",trial,"trials."))
  }
  
  #shift back longitudes to the reference of the original longitudinal range 
  R4coord[,1] = (R4coord[,1] + Lon.range.orig[1] + 360)%%360
  R4coord[R4coord[,1] > 180, 1] = R4coord[R4coord[,1] > 180, 1] - 360
  
  if(trial == MAXtrial) R4coord = NULL
  R4coord
  #list(Rcoord=R4coord, angles=c(angle1, angle2, angle3, mir))
}

##########################################################################################################################
# rrotGlobe() performs multiple random rotations of latitude and longitude coordinates around at least three random axes 
# of a sphere, possibly until the set of coordinates fall within a limited window. Adapted to situations where the range 
# of latitude and longitude is not too small (large or unlimited window), and performs well even when the set of coordinates 
# to move are widely distributed around the globe (i.e. spanning more than 90° longitude and latitude). 

rrotGlobe = function(coord, Lon.range=c(-180,180), Lat.range=c(-90,90), mirror=c("random"), 
                     verbose=T, MAXtrial=1000, MINtrial=3) 
{
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  Rcoord = coord
  
  #if requested, take mirror image along longitude and recenter on original coordinates
  mir=F
  if(mirror==T | substr(mirror,1,1)=="y" | (substr(mirror,1,1)=="r" & runif(1)<0.5)) mir=T 
  if(mir) Rcoord[,1] = -Rcoord[,1] + max(Rcoord[,1]) + min(Rcoord[,1])
  
  #operate of longitudinal shift so that Lon.range starts at 0 
  Lon.range.orig = Lon.range
  Rcoord[,1] = (Rcoord[,1] - Lon.range.orig[1] + 360)%%360 
  Lon.range = (Lon.range - Lon.range.orig[1] + 360)%%360 
  if(Lon.range[2] == 0) Lon.range[2] = 360
  #operate random rotations until points fall within the lat/lon ranges
  for(trial in 1:MAXtrial){
    #chose a random axis (uniform distribution around a sphere)
    axis = c((-acos(runif(1, min=-1, max=1))*180/pi + 90), runif(1)*360)
    Rcoord = RotLonLatAxis(Rcoord, axis, runif(1)*360)
    Rcoord[,1] = (Rcoord[,1] + 360)%%360 
    if(trial>=MINtrial & min(Rcoord[,2])>=Lat.range[1] & max(Rcoord[,2])<=Lat.range[2] & 
        min(Rcoord[,1])>=Lon.range[1] & max(Rcoord[,1])<=Lon.range[2] ) break
  }
  
  if(verbose){
    if( min(Rcoord[,2]) < Lat.range[1] | max(Rcoord[,2]) > Lat.range[2] |  
        min(Rcoord[,1]) < Lon.range[1] | max(Rcoord[,1]) > Lon.range[2]   )
        print(paste("Warning: failed to move original coordinates into the latitudinal or longitudinal ranges after",trial,"trials."))
  }
  
  
  
  #shift back longitudes to the reference of the original longitudinal range 
  Rcoord[,1] = (Rcoord[,1] + Lon.range.orig[1] + 360)%%360 
  Rcoord[Rcoord[,1] > 180, 1] = Rcoord[Rcoord[,1] > 180, 1] - 360
  
  if(trial == MAXtrial) Rcoord = NULL
  Rcoord
  #list(Rcoord=Rcoord, Ntrial=i)
}




