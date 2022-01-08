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
  obj   <- list(x=1:xin, y=1:yin, z = z)
  tempx <- seq(1,xin,length.out=xout)
  tempy <- seq(1,yin,length.out=yout)
  loc   <- make.surface.grid(list(tempx,tempy))
  look  <- interp.surface(obj,loc)
  return(as.surface(loc,look)$z)
}


##--LOAD DATA--#################################
load_data <- function(modis,oc,bathy){
  lats <<- seq(-90,90,length.out=720)
  lons <<- seq(-180,180,length.out=1440)
  
  latsi <<- which(lats > 35 & lats < 80)
  lonsi <<- which(lons > -75 & lons < 35)
  nlat  <<- length(latsi)
  nlon  <<- length(lonsi)
  
  if(modis==TRUE){
    load('data/chl.rdata',envir=.GlobalEnv)
    load('data/sst.rdata',envir=.GlobalEnv)
    load('data/time.rdata',envir=.GlobalEnv)
    load('data/date.rdata',envir=.GlobalEnv)
  }
  if(oc==TRUE){
    load('data/chl_oc.rdata')
    load('data/time_oc.rdata')
    load('data/date_oc.rdata')
    chl_oc  <<- chl_oc[,,224:1070] 
    date_oc <<- date_oc[224:1070] 
    time_oc <<- time_oc[224:1070]
  }
  if(bathy==TRUE){
    file  <- 'data/GEBCO_BATHY_2002-01-01_rgb_1440x720.csv' 
    bathy <- as.matrix(read.csv(file,header=FALSE))
    bathy[bathy==99999] <- NA
    bathy     <- -t(bathy)[,720:1]
    bathy_atl <<- bathy[lonsi,latsi]
  }
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


###############################################################################
## ANALYSIS FUNCTION ##########################################################
###############################################################################
fit <- function(chl,nlon,nlat,date,time,oc){
  location     =shape     =scale     =location_se     =shape_se     =scale_se     =record=record_se=cormap       <- matrix(NA,nrow=nlon,ncol=nlat)
  location_half=shape_half=scale_half=location_half_se=shape_half_se=scale_half_se <- matrix(NA,nrow=nlon,ncol=nlat)
  location_doub=shape_doub=scale_doub=location_doub_se=shape_doub_se=scale_doub_se <- matrix(NA,nrow=nlon,ncol=nlat)
  
  location_trnd=scale_trnd=shape_trnd           <- matrix(NA,nrow=nlon,ncol=nlat)
  location_trnd_se=scale_trnd_se=shape_trnd_se  <- matrix(NA,nrow=nlon,ncol=nlat)
  
  fitBIC=fitlocBIC=fitsclBIC=fitshpBIC          <- matrix(NA,nrow=nlon,ncol=nlat)
  
  maxmonth <- matrix(NA,nrow=nlon,ncol=nlat)
  
  chl_trend=sst_trend <- matrix(NA,nrow=nlon,ncol=nlat)
  
  N_half    <- round(length(date)/2)
  date_half <- date[1:N_half]
  date_doub <- seq.Date(ymd(date[1]), ymd("2040-11-25"),8)
  N_doub    <- length(date_doub)
  
  periods              <- seq(2,100)
  return               <- array(NA,dim=c(nlon,nlat,length(periods)))
  qq                   <- array(NA,dim=c(nlon,nlat,2,19))
  years                <- 1:19
  
  for(i in 1:nlon){ print(i)
    for(j in 1:nlat){
      chl_tmp <- chl[i,j,]
      sst_tmp <- sst[i,j,]
      
      if(sum(!is.na(chl_tmp))>200){
        chl_tmp[is.na(chl_tmp)] <- 0
        
        #compute maxima and month
        proc          <- block_maxima(chl_tmp,date=date)
        maxima        <- proc$maxima
        month         <- proc$month
        maxima_half   <- maxima[1:10]
        maxima_doub   <- c(maxima,maxima)
        
        #find most frequent bloom month
        maxmonth[i,j] <- Mode(month)[1]      
        
        #fit stationary and non-stationary GEVDs
        fit           <- fevd(maxima)
        fitloc        <- fevd(maxima,location.fun=~years)
        fitscl        <- fevd(maxima,scale.fun=   ~years)
        fitshp        <- fevd(maxima,shape.fun=   ~years)
        
        #compute standard errors and store the results
        ses <- sqrt(diag(parcov.fevd(fit))) #extract standard errors
        location_se[i,j] <- ses[1]
        shape_se[i,j]    <- ses[2]
        scale_se[i,j]    <- ses[3]
        
        #store location, shape, and scale results
        location[i,j] <- fit$results$par[1] 
        shape[i,j]    <- fit$results$par[2]
        scale[i,j]    <- fit$results$par[3]
        
        if(oc==FALSE){
          #fit stationary GEVD to halved and doubled time series of maxima
          fit_half <- fevd(maxima_half)
          fit_doub <- fevd(maxima_doub)
          
          location_half[i,j] <- fit_half$results$par[1] 
          shape_half[i,j]    <- fit_half$results$par[2]
          scale_half[i,j]    <- fit_half$results$par[3]
          
          location_doub[i,j] <- fit_doub$results$par[1] 
          shape_doub[i,j]    <- fit_doub$results$par[2]
          scale_doub[i,j]    <- fit_doub$results$par[3]
          
          #store trend results
          location_trnd[i,j]    <- fitloc$results$par[2]
          loc_trnd_try          <- try(sqrt(diag(parcov.fevd(fitloc)))[2],silent=TRUE)
          location_trnd_se[i,j] <- ifelse(class(loc_trnd_try)=='try-error',NA,loc_trnd_try)
          
          scale_trnd[i,j]       <- fitscl$results$par[3]
          scl_trnd_try          <- try(sqrt(diag(parcov.fevd(fitscl)))[3],silent=TRUE)
          scale_trnd_se[i,j]    <- ifelse(class(scl_trnd_try)=='try-error',NA,scl_trnd_try)
          
          shape_trnd[i,j]       <- fitshp$results$par[4]
          shp_trnd_try          <- try(sqrt(diag(parcov.fevd(fitshp)))[4],silent=TRUE)
          shape_trnd_se[i,j]    <- ifelse(class(shp_trnd_try)=='try-error',NA,shp_trnd_try)   
          
          ses_half <- sqrt(diag(parcov.fevd(fit_half))) #extract standard errors
          location_half_se[i,j] <- ses_half[1]
          shape_half_se[i,j]    <- ses_half[2]
          scale_half_se[i,j]    <- ses_half[3]
          
          ses_doub <- sqrt(diag(parcov.fevd(fit_doub))) #extract standard errors
          location_doub_se[i,j] <- ses_doub[1]
          shape_doub_se[i,j]    <- ses_doub[2]
          scale_doub_se[i,j]    <- ses_doub[3]
          
          #compute BIC    
          fitBIC[i,j]    <- 2*fit$results$value    + log(19)*3
          fitlocBIC[i,j] <- 2*fitloc$results$value + log(19)*4
          fitsclBIC[i,j] <- 2*fitscl$results$value + log(19)*4
          fitshpBIC[i,j] <- 2*fitshp$results$value + log(19)*4
        }        
        ##--STORE QUANTILES--##################
        vv <- plot(fit,type='qq')
        
        qq[i,j,1,]  <- vv[,1]  #observed quantiles
        qq[i,j,2,]  <- vv[,2]  #fitted quantiles
        cormap[i,j] <- cor(vv[,1],vv[,2])
        
        if(oc==FALSE){
          if(sum(!is.na(chl_tmp))>200 & sum(!is.na(sst_tmp))>200){
            chl_lm  <- lm(chl_tmp ~ time)    
            sst_lm  <- lm(sst_tmp ~ time)    
            chl_trend[i,j] <- summary(chl_lm)$coefficients[2,1]    
            sst_trend[i,j] <- summary(sst_lm)$coefficients[2,1] 
          }
        }        
      }
    }
  }
  return(list(location=location,           scale=scale,           shape=shape,   #uncomment to organize results from MODIS data
              location_se=location_se,     scale_se=scale_se,     shape_se=shape_se,
              location_trnd=location_trnd, scale_trnd=scale_trnd, shape_trnd=shape_trnd,
              fitlocBIC=fitlocBIC,         fitsclBIC=fitsclBIC,   fitshpBIC=fitshpBIC,   fitBIC=fitBIC,
              record=record, 
              record_se=record_se,
              cormap=cormap,
              chl_trend=chl_trend, sst_trend=sst_trend,
              qq=qq,
              location_doub=location_doub,  scale_doub=scale_doub,  shape_doub=shape_doub,
              location_half=location_half,  scale_half=scale_half,  shape_half=shape_half,
              location_doub_se=location_doub_se,  scale_doub_se=scale_doub_se,  shape_doub_se=shape_doub_se,
              location_half_se=location_half_se,  scale_half_se=scale_half_se,  shape_half_se=shape_half_se,
              maxmonth=maxmonth))
}
