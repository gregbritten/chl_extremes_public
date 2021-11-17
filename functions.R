##--PLOTTING FUNCTION--#####################
image.plot2 <- function(x,y,z,zlim,cols){
  tmp <- z
  tmp[tmp<zlim[1]] <- zlim[1]
  tmp[tmp>zlim[2]] <- zlim[2]
  image.plot(x=x,y=y,z=tmp,zlim=zlim,col=cols,xaxt='n',yaxt='n',ylab='',xlab='')
}


##--BINLINEAR INTERPOLATION--#####################
resize_bilinear <- function(xin,yin,xout,yout,z){
  library(fields)
  obj <- list(x=1:xin, y=1:yin, z = z)
  tempx <- seq(1,xin,xout)
  tempy <- seq(1,yin,yout)
  loc <- make.surface.grid(list(tempx,tempy))
  
  look <- interp.surface( obj, loc)
  return(as.surface( loc, look)$z)
}


##--FUNCTION TO EXTRACT BLOCK MAXIMA--###################################
block_maxima <- function(x,date){
  years <- year(date)
  months <- month(date)
  year_ind <- unique(years)
  maxima=month <- numeric(length(year_ind))
  for(i in 1:length(year_ind)){
    mmonths   <- months[years==year_ind[i]]
    maxima[i] <- max(x[years==year_ind[i]],na.rm=TRUE)
    month[i]  <- mmonths[x[years==year_ind[i]]==maxima[i]][1]
  }
  return(list(maxima=maxima,month=month))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
