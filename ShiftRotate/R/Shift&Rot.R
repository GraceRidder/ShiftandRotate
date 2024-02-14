LatLon2XYZ = function(coord){
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  coord = coord/180
  a1 = sinpi(1/2 - coord[, 1])
  XYZ = cbind(a1 * cospi(coord[, 2]), a1 * sinpi(coord[, 2]), cospi(1/2-coord[, 1]))
  colnames(XYZ) = c("X", "Y", "Z")
  XYZ
}

XYZ2LatLon = function(XYZ){
  XYZ = as.matrix(XYZ)
  if (ncol(XYZ) == 1) XYZ = t(XYZ)
  Lat = 90 - (acos(XYZ[,3]) * 180/pi)
  Lon = ( (atan(XYZ[, 2]/XYZ[, 1]) + pi*I(XYZ[,1] < 0))%%(2*pi)) * 180/pi
  Lon[XYZ[,3]==1] = 0
  Lon[Lon>180] = Lon[Lon>180] - 360
  cbind(Lat,Lon)
}


RotLatLongAxis = function(coord, axis, angle){
  R_n_a = function(a, n1, n2, n3)  rbind( c( cos(a)+n1^2*(1-cos(a)), n1*n2*(1-cos(a))-n3*sin(a), n1*n3*(1-cos(a))+n2*sin(a)),
                                          c( n1*n2*(1-cos(a))+n3*sin(a), cos(a)+n2^2*(1-cos(a)), n2*n3*(1-cos(a))-n1*sin(a)),
                                          c( n1*n3*(1-cos(a))-n2*sin(a), n2*n3*(1-cos(a))+n1*sin(a), cos(a)+n3^2*(1-cos(a)))  )
  XYZ = LatLon2XYZ(coord)
  xyz = LatLon2XYZ(axis)
  Rna = R_n_a(-1*angle/180*pi, xyz[1], xyz[2], xyz[3])
  RXYZ = XYZ %*% Rna
  Rcoord = XYZ2LatLon(RXYZ)
  Rcoord
}



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



shiftrotGlobe = function(coord, Lat.range=c(-90,90), Lon.range=c(-180,180), rotation=T, Lat.shift=T, Lon.shift=T,
                         mirror=c("random"), verbose=T, MAXtrial=1000)
{
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  Lat = coord[,1]
  Lon = coord[,2]
  Lon.range.orig = Lon.range
  #operate longitudinal shift so that Lon.range starts at 0
  Lon = (Lon - Lon.range.orig[1] + 360)%%360
  Lon.range = (Lon.range - Lon.range[1] + 360)%%360
  if(Lon.range[2] == 0) Lon.range[2] = 360
  #check if points are within lat-long ranges (this does not impede the rest of the algorithm).
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

  Rcoord = cbind(RLat, RLon)

  trial = 0
  repeat{
    trial = trial+1

    #Step 2: apply rotation around axis pointing to centroid of data points, checking that points can enter the window size
    #axis1 = centroid axis (mean vector of coordinates, scaled to unit length)
    axis1 = XYZ2LatLon( colMeans(LatLon2XYZ(Rcoord)) / sqrt(sum(colMeans(LatLon2XYZ(Rcoord))^2)) )
    if(rotation){
      if(!is.logical(rotation) & is.numeric(rotation)){
        R2coord = RotLatLongAxis(Rcoord, axis1, rotation)
      }else{
        angle1 = runif(1)*360
        R2coord = RotLatLongAxis(Rcoord, axis1, angle1)
      }
    } else R2coord = Rcoord

    #Step 3: apply latitudinal rotation along an equatorial vector perpendicular to the centroid vector
    if(Lat.shift){
      axis2 = c(0, (axis1[,2] - 90 + 360)%%360 ) #equatorial axis (lat==0) pointing perpendicular to centroid axis1
      if(!is.logical(Lat.shift) & is.numeric(Lat.shift)){
        R3coord = RotLatLongAxis(R2coord, axis2, Lat.shift)
      }else{
        angle2 = -acos(runif(1, min=sin(Lat.range[1]/180*pi), max=sin(Lat.range[2]/180*pi)))*180/pi + 90 #random latitude for a point uniformly distributed on a sphere between 2 latitudes
        angle2 = angle2 - axis1[1] #rescale angle so that a change in latitude from the cebtroid latitude leads to a probability of final latitude ~ cos(latitude)
        R3coord = RotLatLongAxis(R2coord, axis2, angle2)
      }
    } else R3coord = R2coord

    #Step 4: apply longitudinal translation, adding a random longitudinal shift
    R4coord = R3coord
    if(Lon.shift){
      if(!is.logical(Lon.shift) & is.numeric(Lon.shift)){
        R4coord[,2] = (R4coord[,2] + Lon.shift + 360)%%360
      }else{
        angle3 = 0
        if((Lon.range[2]-Lon.range[1])==360) angle3 = runif(1, min = 0, max=360)
        if((Lon.range[2]-Lon.range[1])<360 & (max(R3coord[,2])-min(R3coord[,2])) <= (Lon.range[2]-Lon.range[1]) )
          angle3 = runif(1, min = Lon.range[1]-min(R3coord[,2]), max = Lon.range[2]-max(R3coord[,2]))
        R4coord[,2] = (R4coord[,2] + angle3 + 360)%%360
      }
    }

    # Finally, check that moved coordinates fall in the window, otherwise try again, unless there is no random moves or MAXtrial is reached
    if(min(R4coord[,1])>=Lat.range[1] & max(R4coord[,1])<=Lat.range[2] &
       min(R4coord[,2])>=Lon.range[1] & max(R4coord[,2])<=Lon.range[2]) break
    if(rotation!=T & Lat.shift!=T & Lon.shift!=T) break
    if(trial == MAXtrial) break
  }

  if(verbose){
    if( min(R4coord[,1]) < Lat.range[1] | max(R4coord[,1]) > Lat.range[2] |
        min(R4coord[,2]) < Lon.range[1] | max(R4coord[,2]) > Lon.range[2]   )
        print(paste("Warning: failed to move original coordinates into the latitudinal or longitudinal ranges after",trial,"trials."))
  }

  #shift back longitudes to the reference of the original longitudinal range
  R4coord[,2] = (R4coord[,2] + Lon.range.orig[1] + 360)%%360
  R4coord[R4coord[,2] > 180, 2] = R4coord[R4coord[,2] > 180, 2] - 360

  if(trial == MAXtrial) NULL
  else R4coord
  #list(Rcoord=R4coord, angles=c(angle1, angle2, angle3, mir))
}


rrotGlobe = function(coord, Lat.range=c(-90,90), Lon.range=c(-180,180), mirror=c("random"),
                     verbose=T, MAXtrial=1000, MINtrial=3)
{
  coord = as.matrix(coord)
  if (ncol(coord) == 1) coord = t(coord)
  Rcoord = coord

  #if requested, take mirror image along longitude and recenter on original coordinates
  mir=F
  if(mirror==T | substr(mirror,1,1)=="y" | (substr(mirror,1,1)=="r" & runif(1)<0.5)) mir=T
  if(mir) Rcoord[,2] = -Rcoord[,2] + max(Rcoord[,2]) + min(Rcoord[,2])

  #operate of longitudinal shift so that Lon.range starts at 0
  Lon.range.orig = Lon.range
  Rcoord[,2] = (Rcoord[,2] - Lon.range.orig[1] + 360)%%360
  Lon.range = (Lon.range - Lon.range.orig[1] + 360)%%360
  if(Lon.range[2] == 0) Lon.range[2] = 360
  #operate random rotations until points fall within the lat/lon ranges
  for(trial in 1:MAXtrial){
    #chose a random axis (uniform distribution around a sphere)
    axis = c((-acos(runif(1, min=-1, max=1))*180/pi + 90), runif(1)*360)
    Rcoord = RotLatLongAxis(Rcoord, axis, runif(1)*360)
    Rcoord[,2] = (Rcoord[,2] + 360)%%360
    if(trial>=MINtrial & min(Rcoord[,1])>=Lat.range[1] & max(Rcoord[,1])<=Lat.range[2] &
        min(Rcoord[,2])>=Lon.range[1] & max(Rcoord[,2])<=Lon.range[2] ) break
  }

  if(verbose){
    if( min(Rcoord[,1]) < Lat.range[1] | max(Rcoord[,1]) > Lat.range[2] |
        min(Rcoord[,2]) < Lon.range[1] | max(Rcoord[,2]) > Lon.range[2]   )
        print(paste("Warning: failed to move original coordinates into the latitudinal or longitudinal ranges after",trial,"trials."))
  }



  #shift back longitudes to the reference of the original longitudinal range
  Rcoord[,2] = (Rcoord[,2] + Lon.range.orig[1] + 360)%%360
  Rcoord[Rcoord[,2] > 180, 2] = Rcoord[Rcoord[,2] > 180, 2] - 360

  if(trial == MAXtrial) Rcoord = NULL
  Rcoord
  #list(Rcoord=Rcoord, Ntrial=i)
}


#test functions: 1 shiftrotXY, 2 shiftrotGlobe, 3 rrotGlobe
fun = 1
library(Directional) #to visualize points on a sphere using sphereplot()

coord=data.frame(lat=c(1,6,11,16,21,26,6),long=c(-80,-80,-80,-80,-80,-80,-75))
Lat.range = c(-20,40)
Lon.range = c(-90,-10)
rotation = Lat.shift = Lon.shift = T
mirror = "r"

plot(coord[,2:1], asp=1, xlim=c(-180,180), ylim=c(-90,90))
lines( cbind(c(Lon.range[1],Lon.range[2],Lon.range[2],Lon.range[1],Lon.range[1]), c(Lat.range[1],Lat.range[1],Lat.range[2],Lat.range[2],Lat.range[1])))
sphereplot(LatLon2XYZ(coord), col=1) #projection on a sphere
for(i in 1:6){
  if(fun==1) newcoord = shiftrotXY(coord[,2:1], Y.range=Lat.range, X.range=Lon.range, Y.shift=Lat.shift, X.shift=Lon.shift, rotation=rotation, mirror=mirror)[,2:1]
  if(fun==2) newcoord = shiftrotGlobe(coord, Lat.range=Lat.range, Lon.range=Lon.range, Lat.shift=Lat.shift, Lon.shift=Lon.shift, rotation=rotation, mirror=mirror)
  if(fun==3) newcoord = rrotGlobe(coord, Lat.range=Lat.range, Lon.range=Lon.range, mirror=mirror)
  points(newcoord[,2:1], col=1+i)
#  sphereplot(LatLon2XYZ(newcoord), col=1+i)
}


#try also with
Lat.range = c(-90,90); Lon.range = c(-180,180) #unlimited lat/lon ranges
Lat.range = c(30,70); Lon.range = c(-90,10)    #window outside original lat/lon ranges
Lat.range = c(-70,0); Lon.range = c(160,-160)  #longitude range crossing the 180°/-180° limit
Lat.range = c(0,10); Lon.range = c(0,10)  #too small window

coord=data.frame(lat=c(1,6,11,16,21,26,6),long=c(178,178,178,178,178,178,-177)) #points across the 180°/-180° limit
coord=data.frame(lat=c(41,46,51,56,61,66,46,41,46,51),long=c(-80,-80,-80,-80,-80,-80,-75,80,80,80)) #points dispersed on opposite side of the world => not ideal for shiftrotGlobe()
Lat.range = range(coord[,1]) + c(-20,20); Lon.range = range(coord[,2]) + c(-20,20)

#and with
rotation = F

#check distribution of lat and long and inter-point distances after randomization
coord=cbind(lat=c(50,50.0001),lon=c(0,0.0001))

rotation = Lat.shift = Lon.shift = T
mirror = "r"
gdist = distm(coord[,2:1], fun=distCosine) #distances btw points

plot(coord[,2:1], asp=1, xlim=c(-180,180), ylim=c(-90,90))
lines( cbind(c(Lon.range[1],Lon.range[2],Lon.range[2],Lon.range[1],Lon.range[1]), c(Lat.range[1],Lat.range[1],Lat.range[2],Lat.range[2],Lat.range[1])))
Nrep=5000
newcoordA = array(dim=c(dim(coord),Nrep))
gdistA = array(dim=c(dim(gdist),Nrep))
for(i in 1:Nrep){
  if(fun==1) R = shiftrotXY(coord[,2:1], Y.range=Lat.range, X.range=Lon.range, Y.shift=Lat.shift, X.shift=Lon.shift, rotation=rotation, mirror=mirror)[,2:1]
  if(fun==2) R = shiftrotGlobe(coord, Lat.range=Lat.range, Lon.range=Lon.range, Lat.shift=Lat.shift, Lon.shift=Lon.shift, rotation=rotation, mirror=mirror)
  if(fun==3) R = rrotGlobe(coord, Lat.range=Lat.range, Lon.range=Lon.range, mirror=mirror)
  newcoordA[,,i] = R
  points(R[,2:1], col=1+i)
}
hist(newcoordA[,1,], nclass=36) #Latitude => we expect less points close to borders or at higher |latitude|
range(newcoordA[,1,])
hist(newcoordA[,2,], 36) #Longitude => we expect less points close to borders or uniform distribution if Lon.range = 360°
range(newcoordA[,2,])
#expected distribution across latitude without constrain
plot(x=(-90:90), y=cos((-90:90)/180*pi))

hist(gdistA)
range(gdistA)

# Rotations by systematic angles
coord=data.frame(lat=c(1,6,11,16,21,26,6),long=c(-80,-80,-80,-80,-80,-80,-75))
Lat.range = c(-90,90)
Lon.range = c(-180,180)
plot(coord[,2:1], asp=1, xlim=Lon.range, ylim=Lat.range)
lines( cbind(c(Lon.range[1],Lon.range[2],Lon.range[2],Lon.range[1],Lon.range[1]), c(Lat.range[1],Lat.range[1],Lat.range[2],Lat.range[2],Lat.range[1])))
rotation = Lat.shift = Lon.shift = mirror = F
for(i in 1:7){
  rotation = i*45
  newcoord = shiftrotGlobe(coord,Lat.range=Lat.range, Lon.range=Lon.range, Lat.shift=Lat.shift, Lon.shift=Lon.shift, rotation=rotation, mirror=mirror)
  points(newcoord[,2:1], col=1+i)
}
coord2 = newcoord
rotation = F
for(i in 1:18){
  Lon.shift = i*20
  newcoord = shiftrotGlobe(coord2,Lat.range=Lat.range, Lon.range=Lon.range, Lat.shift=Lat.shift, Lon.shift=Lon.shift, rotation=rotation, mirror=mirror)
  points(newcoord[,2:1], col=1+i)
}
Lon.shift = F
for(i in 1:18){
  Lat.shift = i*20
  newcoord = shiftrotGlobe(coord2,Lat.range=Lat.range, Lon.range=Lon.range, Lat.shift=Lat.shift, Lon.shift=Lon.shift, rotation=rotation, mirror=mirror)
  points(newcoord[,2:1], col=1+i)
}


# illustration of rotation+shifts step by step
coord=data.frame(lat=c(1,2,3,4,5,6,2),long=c(1,1,1,1,1,1,1.5))

plot(coord[,2:1], asp=1, xlim=c(min(coord[,2])-10,max(coord[,2])+10), ylim=c(min(coord[,1])-10,max(coord[,1])+10))
sphereplot(LatLon2XYZ(coord))

#rotate around centroid
angle=30
axis=colMeans(coord)
Rcoord=RotLatLongAxis(coord, axis, angle)
points(Rcoord[,2:1], col=2)

#rotate along polar axis (longitude)
angle=5
axis=c(90,0)
R2coord=RotLatLongAxis(Rcoord, axis, angle)
points(R2coord[,2:1], col=3)

#rotate along equatorial axis perpendicular to centroid
angle=5
axis=c(0,mean(coord[,2])-90)
R3coord=RotLatLongAxis(R2coord, axis, angle)
points(R3coord[,2:1], col=4)


