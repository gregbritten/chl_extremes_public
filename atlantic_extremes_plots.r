rm(list=ls())

library(maps)
setwd('~/dropbox/working/chlorophyll/data/')

load('ATL_DAT.rdata')
DAT <- DAT
# load('~/dropbox/working/chlorophyll/data/ATL_DAT_oc.rdata')
# DAT <- DAT_oc

lats <- seq(-90,90,length.out=720)
lons <- seq(-180,180,length.out=1440)

latsi <- which(lats > 35 & lats < 80)
lonsi <- which(lons > -75 & lons < 35)
nlat  <- length(latsi)
nlon  <- length(lonsi)

file  <- 'GEBCO_BATHY_2002-01-01_rgb_1440x720.csv' 
bathy <- as.matrix(read.csv(file,header=FALSE))
bathy[bathy==99999] <- NA
bathy     <- -t(bathy)[,720:1]
bathy_atl <- bathy[lonsi,latsi]

##--PLOTTING FUNCTION--#####################
image.plot2 <- function(x,y,z,zlim,cols){
  tmp <- z
  tmp[tmp<zlim[1]] <- zlim[1]
  tmp[tmp>zlim[2]] <- zlim[2]
  image.plot(x=x,y=y,z=tmp,zlim=zlim,col=cols,xaxt='n',yaxt='n',ylab='',xlab='')
}

# Xloc <- data.frame(mod=c(DAT$location),occi=c(DAT_oc$location))
# Xscl <- data.frame(mod=c(DAT$scale),   occi=c(DAT_oc$scale))
# Xshp <- data.frame(mod=c(DAT$shape),   occi=c(DAT_oc$shape))
# 
# #remove outliers from pathological fits
# Xloc$mod[Xloc$mod<0]   = Xloc$occi[Xloc$occi<0]   = Xloc$mod[Xloc$mod>10] = Xloc$occi[Xloc$occi>10] <- NA
# Xscl$mod[Xscl$mod< -2] = Xscl$occi[Xscl$occi< -2] = Xscl$mod[Xscl$mod>4]  = Xscl$occi[Xscl$occi>4] <- NA
# Xshp$mod[Xshp$mod<0]   = Xshp$occi[Xshp$occi<0]   = Xshp$mod[Xshp$mod>8]  = Xshp$occi[Xshp$occi>8] <- NA
# 
# Xloc <- Xloc[complete.cases(Xloc),]
# Xscl <- Xscl[complete.cases(Xscl),]
# Xshp <- Xshp[complete.cases(Xshp),]
# 
# samp <- sample(size=500,1:nrow(Xloc))  
# par(mfrow=c(1,3))
# plot(Xloc$mod[samp],Xloc$occi[samp],ylim=c(0,4),xlim=c(0,4))
# plot(Xscl$mod[samp],Xscl$occi[samp],ylim=c(-1,2),xlim=c(-1,2))
# plot(Xshp$mod[samp],Xshp$occi[samp],ylim=c(0,4),xlim=c(0,4))

##########################################################################  
##--MONTH OF BLOOM--######################################################
##########################################################################  
pdf('~/dropbox/working/chlorophyll/plots/maxmonth.pdf',height=6,width=7)
par(mfrow=c(1,1), oma=c(2,2,2,2))
image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$maxmonth,zlim=c(1,12),col=turbo(12))
map(add=TRUE,col='grey',fill=TRUE)
axis(side=1); axis(side=2)
box(lwd=3)
mtext(expression('Bloom Month'),adj=0,line=0.5)
mtext(side=1,line=2.5,'Longitude')
mtext(side=2,line=2.5,'Latitude')
dev.off()



##########################################################################  
##--MAP FITTED PARAMETERS--###############################################
##########################################################################  
cols <- turbo(20)

#pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_maps_revision.pdf',height=4.5,width=14)
pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_maps_oc_revision.pdf',height=4.5,width=14)
par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))
  
  # image.plot2(x=lons[lonsi],y=lats[latsi],z=bathy_atl/1000,zlim=c(0,5),col=cols)
  #   axis(side=1,labels=NA); axis(side=2)
  #   #contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
  #   map(add=TRUE,col='grey',fill=TRUE)
  #   box(lwd=3)
  #   mtext('a) Bathymetry',adj=0,line=0.5)
  #   mtext('[km]',adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location,zlim=c(0,4),col=cols)
    contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2)
    box(lwd=3)
    mtext(expression('a) Location'~mu),adj=0,line=0.5)
    mtext(expression('[mgChl/m'^3*']'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale,zlim=c(-0.5,2),col=cols)
    contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2,labels=NA)
    box(lwd=3)
    mtext(expression('b) Scale'~sigma),adj=0,line=0.5)
    mtext(expression('[mgChl/m'^3*']'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape,zlim=c(0,2.5),col=cols)
    contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2,labels=NA)
    box(lwd=3)
    axis(side=1); axis(side=2,labels=NA)
    mtext(expression('c) Shape'~xi),adj=0,line=0.5)
    mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  mtext('Latitude',outer=TRUE,side=2,line=1.5,cex=1.2)
  mtext('Longitude',outer=TRUE,side=1,line=1.5,cex=1.2)
dev.off()

##########################################################################  
##--FITTED PARAMETER SES--################################################
########################################################################## 
cols <- turbo(20)

#pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_se_maps_oc_revision.pdf',height=4.5,width=14)
pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_se_maps_revision.pdf',height=4.5,width=14)
par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))

  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_se/DAT$location,zlim=c(0,0.5),col=cols)
    #contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2)
    box(lwd=3)
    mtext(expression('a) sd('*mu*')/'*mu),adj=0,line=0.5)
    mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_se/DAT$scale,zlim=c(0,2),col=cols)
    #contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2,labels=NA)
    box(lwd=3)
    mtext(expression('b) sd('*sigma*')/'*sigma),adj=0,line=0.5)
    mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_se/DAT$shape,zlim=c(0,0.5),col=cols)
    #contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2,labels=NA)
    box(lwd=3)
    axis(side=1); axis(side=2,labels=NA)
    mtext(expression('d) sd('*xi*')/'*xi),adj=0,line=0.5)
    mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  mtext('Latitude',outer=TRUE,side=2,line=1.5,cex=1.2)
  mtext('Longitude',outer=TRUE,side=1,line=1.5,cex=1.2)
dev.off()



##########################################################################  
##--POWER ANALYSIS--######################################################
########################################################################## 
cols <- turbo(20)

pdf('~/dropbox/working/chlorophyll/plots/gevd_parameter_se_power.pdf',height=4*2,width=12)
par(mfrow=c(3,3),mar=c(1,1,2,4),oma=c(4,4,2,2))

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_se/DAT$location,zlim=c(0,0.5),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('a) sd('*mu*')/'*mu),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_se/DAT$scale,zlim=c(0,2),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('b) sd('*sigma*')/'*sigma),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_se/DAT$shape,zlim=c(0,0.5),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  axis(side=1); axis(side=2,labels=NA)
  mtext(expression('c) sd('*xi*')/'*xi),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_half_se/DAT$location_half,zlim=c(0,0.5),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2)
  box(lwd=3)
  mtext(expression('d)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_half_se/DAT$scale_half,zlim=c(0,2),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  #mtext('Halving Sample Size',line=0.5)
  mtext(expression('e)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_half_se/DAT$shape_half,zlim=c(0,0.5),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('f)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_doub_se/DAT$location_doub,zlim=c(0,0.5),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('g)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_doub_se/DAT$scale_doub,zlim=c(0,2),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  #mtext('Doubling Sample Size',line=0.5)
  mtext(expression('h)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_doub_se/DAT$shape_doub,zlim=c(0,0.5),col=cols)
  map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('i)'),adj=0,line=0.5)

    mtext('Latitude',outer=TRUE,side=2,line=1.5,cex=1.2)
    mtext('Longitude',outer=TRUE,side=1,line=1.5,cex=1.2)
dev.off()






################################################################
##--BINNED ANALYSIS--###########################################
################################################################
##--BATHYMETRY SCATTER--#####################
pdf('~/dropbox/working/chlorophyll/plots/parameter_bathymetry_revision.pdf',height=6,width=9)
par(mfrow=c(2,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
  location2 <- DAT$location
  location2[location2>5] <- 5
  location_bins <- stats.bin(x=c(bathy_atl),y=c(location2),breaks=seq(50,5000,length.out=500))
  location_xx <- cbind(location_bins$centers,location_bins$stats[2,],location_bins$stats[3,],location_bins$stats[1,])
  location_xx <- location_xx[complete.cases(location_xx),]
  plot(location_xx[,1],location_xx[,2],type='l',lwd=2)
  lines(location_xx[,1],location_xx[,2] + 2*location_xx[,3]/sqrt(location_xx[,4]),type='l',lty=2)#,col='grey')
  lines(location_xx[,1],location_xx[,2] - 2*location_xx[,3]/sqrt(location_xx[,4]),type='l',lty=2)#,col='grey')
  #lines(location_xx[,1],location_xx[,2] + location_xx[,3],type='l',lty=2,col='red')
  #lines(location_xx[,1],location_xx[,2] - location_xx[,3],type='l',lty=2,col='red')
  mtext(side=2,expression('Location'~mu~'[mgChl/m'^3*']'),line=2.25,cex=0.8)
  mtext(side=1,'Bathymetry [m]',line=2.5,cex=0.8)
  mtext('a)',adj=-0.1)
  
  scale2 <- DAT$scale
  scale2[scale2>5] <- 5
  scale_bins <- stats.bin(x=c(bathy_atl),y=c(scale2),breaks=seq(50,5000,length.out=500))
  scale_xx <- cbind(scale_bins$centers,scale_bins$stats[2,],scale_bins$stats[3,],scale_bins$stats[1,])
  scale_xx <- scale_xx[complete.cases(scale_xx),]
  plot(scale_xx[,1],scale_xx[,2],type='l',lwd=2)
  lines(scale_xx[,1],scale_xx[,2] + 2*scale_xx[,3]/sqrt(scale_xx[,4]),type='l',lty=2)#,col='grey')
  lines(scale_xx[,1],scale_xx[,2] - 2*scale_xx[,3]/sqrt(scale_xx[,4]),type='l',lty=2)#,col='grey')
  mtext(side=2,expression('Scale'~sigma~'[mgChl/m'^3*']'),line=2.25,cex=0.8)
  mtext(side=1,'Bathymetry [m]',line=2.5,cex=0.8)
  mtext('b)',adj=-0.1)
  
  shape2 <- DAT$shape
  shape2[shape2>5] <- 5
  shape_bins <- stats.bin(x=c(bathy_atl),y=c(shape2),breaks=seq(50,5000,length.out=500))
  shape_xx <- cbind(shape_bins$centers,shape_bins$stats[2,],shape_bins$stats[3,],shape_bins$stats[1,])
  shape_xx <- shape_xx[complete.cases(shape_xx),]
  plot(shape_xx[,1],shape_xx[,2],type='l',lwd=2)
  lines(shape_xx[,1],shape_xx[,2] + 2*shape_xx[,3]/sqrt(shape_xx[,4]),type='l',lty=2)#,col='grey')
  lines(shape_xx[,1],shape_xx[,2] - 2*shape_xx[,3]/sqrt(shape_xx[,4]),type='l',lty=2)#,col='grey')
  mtext(side=2,'Shape'~xi~'[unitless]',line=2.25,cex=0.8)
  mtext(side=1,'Bathymetry [m]',line=2.5,cex=0.8)
  mtext('c)',adj=-0.1)
  
  tmpxx <- data.frame(x=location_xx[,2],y=scale_xx[,2])
  tmpxx <- tmpxx[tmpxx[,1]<2,]
  plot(tmpxx[,1],tmpxx[,2],pch=16,cex=0.7)
  fit <- lm(y ~ x, data=tmpxx)
  lines(seq(0.6,1.8,0.01),predict(fit, newdata=list(x=seq(0.6,1.8,0.01))))
  mtext(side=1,expression('Location'~mu~'[mgChl/m'^3*']'),line=2.5,cex=0.8); 
  mtext(side=2,expression('Scale'~sigma~'[mgChl/m'^3*']'),line=2.25,cex=0.8)
  mtext('d)',adj=-0.1,line=0.1)
  int <- round(summary(fit)$coefficients[1,1],2)
  slp <- round(summary(fit)$coefficients[2,1],2)
  legend('topleft',legend=bquote(italic('y')*' = '*.(slp)*italic('x')*' + '*.(int)),bty='n')
  
  
  tmpxx <- data.frame(x=location_xx[,2],y=shape_xx[,2])
  tmpxx <- tmpxx[tmpxx[,1]<2,]
  plot(tmpxx[,1],tmpxx[,2],pch=16,cex=0.7)
  fit <- lm(y ~ x, data=tmpxx)
  lines(seq(0.65,1.9,0.01),predict(fit, newdata=list(x=seq(0.65,1.9,0.01))))
  mtext(side=1,expression('Location'~mu~'[mgChl/m'^3*']'),line=2.5,cex=0.8); 
  mtext(side=2,expression('Shape'~xi~'[unitless]'),line=2.25,cex=0.8)
  mtext('e)',adj=-0.1,line=0.1)
  int <- round(summary(fit)$coefficients[1,1],2)
  slp <- round(summary(fit)$coefficients[2,1],2)
  legend('topleft',legend=bquote(italic('y')*' = '*.(slp)*italic('x')*' + '*.(int)),bty='n')
  
  tmpxx <- data.frame(x=shape_xx[,2],y=scale_xx[,2])
  tmpxx <- tmpxx[tmpxx[,1]<1,]
  plot(tmpxx[,1],tmpxx[,2],pch=16,cex=0.7)
  fit <- lm(y ~ x, data=tmpxx)
  lines(seq(0.25,0.9,0.01),predict(fit, newdata=list(x=seq(0.25,0.9,0.01))))
  mtext(side=1,expression('Shape'~xi~'[unitess]'),line=2.5,cex=0.8); 
  mtext(side=2,expression('Scale'~sigma~'[mgChl/m'^3*']'),line=2.25,cex=0.8)
  mtext('f)',adj=-0.1,line=0.1)
  int <- round(summary(fit)$coefficients[1,1],2)
  slp <- round(summary(fit)$coefficients[2,1],2)
  legend('topleft',legend=bquote(italic('y')*' = '*.(slp)*italic('x')*' + '*.(int)),bty='n')
dev.off()


#########################################################################################
##--QQ PLOT GOONESS OF FIT--#############################################################
#########################################################################################
X   <- cbind(c(qq[,,1,]),c(qq[,,2,]))
X2  <- X[X[,1]<10 & X[,1]>1E-5 & X[,2]<10 & X[,2]>1E-5,]  #remote extreme outliers where the parameter optimization failed
XX  <- X2[sample(size=5000,1:nrow(X2)),]             #take a random sample for plotting purposes

cor(X2[,1],X2[,2],use='pairwise.complete.obs') #evaluate correlation
cor(log10(X2[,1]),log10(X2[,2]),use='pairwise.complete.obs')

#pdf('~/dropbox/working/chlorophyll/plots/qqplot.pdf',height=4.5, width=9)
pdf('~/dropbox/working/chlorophyll/plots/qqplot_oc.pdf',height=4.5, width=9)
#layout(matrix(c(1,2,3,3),byrow=TRUE,ncol=2))
par(mfrow=c(1,2),mar=c(2,3,2,2),oma=c(2,2,2,2),cex.axis=0.8)
  plot(XX[,1],XX[,2],xlim=c(0,10),ylim=c(0,10),pch=16,cex=0.2,bty='n')
  abline(0,1)
  mtext(side=1,expression('Observed Bloom Quantiles [mgChl/m'^3*']'),line=2.5)
  mtext(side=2,expression('Fitted GEVD Quantiles [mgChl/m'^3*']'),line=2.5)
  mtext('a)',adj=-0.1)
  
  plot(log10(XX[,1]),log10(XX[,2]),xlim=c(-1,2),ylim=c(-1,2),pch=16,cex=0.2,bty='n')
  abline(0,1)
  mtext(side=1,expression('log'['10']~'(Observed Bloom Quantiles [mgChl/m'^3*'])'),line=2.5)
  mtext(side=2,expression('log'['10']~'(Fitted GEVD Quantiles [mgChl/m'^3*'])'),line=2.5)
  mtext('b)',adj=-0.1)
dev.off()

#########################################################################################
##--QQ MAPS--############################################################################
#########################################################################################
#pdf('~/dropbox/working/chlorophyll/plots/qqmap.pdf',height=5,width=9)
pdf('~/dropbox/working/chlorophyll/plots/qqmap_oc.pdf',height=5,width=9)
par(mfrow=c(1,1))
  image.plot(x=lons[lonsi],y=lats[latsi],cormap,col=turbo(20),zlim=c(0.7,1.0),xlab='',ylab='')
  map(add=TRUE)
  box(lwd=2)
  mtext(side=1,'Longitude',line=2.5)
  mtext(side=2,'Latitude',line=2.5)
dev.off()
#########################################################################################
##--NON-STATIONARY & TRENDS--############################################################
#########################################################################################
redblue <- colorRampPalette(c('dark blue','blue','white','red','dark red'))

dBIC_loc <- DAT$fitBIC - DAT$fitlocBIC
dBIC_loc[dBIC_loc<0] <- 0
dBIC_scl <- DAT$fitBIC - DAT$fitsclBIC
dBIC_scl[dBIC_scl<0] <- 0
dBIC_shp <- DAT$fitBIC - DAT$fitshpBIC
dBIC_shp[dBIC_shp<0] <- 0

pdf('~/dropbox/working/chlorophyll/plots/nonstationary_BIC_chl_sst.pdf',height=9,width=10)
par(mfrow=c(3,2),mar=c(2,2,2,3),oma=c(2,2,2,2))
image.plot(x=lons[lonsi],y=lats[latsi],dBIC_loc,zlim=c(0,10),col=turbo(20))
  map(add=TRUE)
  mtext(expression('a) '*Delta*'BIC Location'~mu~'[mgChl/m'^3*']'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],chl_trend,zlim=c(-0.1,0.1),col=redblue(20))
  map(add=TRUE)
  mtext(expression('d) Chlorophyll trend [mgChl/m'^3*'/yr]'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],dBIC_scl,zlim=c(0,10),col=turbo(20))
  map(add=TRUE)
  mtext(expression('b) '*Delta*'BIC Scale'~sigma~'[mgChl/m'^3*']'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],sst_trend,zlim=c(-0.2,0.2),col=redblue(20))
  map(add=TRUE)
  mtext(expression('e) Sea Surface Temperature trend ['*degree*'C/yr]'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],dBIC_shp,zlim=c(0,10),col=turbo(20))
  map(add=TRUE)
  mtext(expression('c) '*Delta*'BIC Shape'~xi~'[unitless]'), adj=0,line=0.25)
  box(lwd=2)
mtext('Latitude',side=2,outer=TRUE,line=0.5)
mtext('Longitude',side=1,outer=TRUE)
dev.off()


##########################################################################
##--BASIN MASKS--#########################################################
##########################################################################
nc <- nc_open('~/dropbox/data/masks/mask.nc')
image.plot(ncvar_get(nc,'atlantic'))


