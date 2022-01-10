rm(list=ls())
gc()

library(maps)
source('functions.r')

load_data(modis=FALSE,oc=FALSE,bathy=TRUE)

##--MAKE PLOTS WITH MODIS--#########################
load('data/ATL_DAT.rdata')
DAT <- DAT
ocyn <- ''

##--MAKE PLOTS WITH OC-CCI--#########################
#load('data/ATL_DAT_oc.rdata')
#DAT <- DAT_oc
#ocyn <- '_oc'

##########################################################################  
##--MONTH OF BLOOM--######################################################
##########################################################################  
pdf(paste0('~/dropbox/working/chlorophyll/plots/maxmonth',ocyn,'.pdf'),height=6,width=7)
par(mfrow=c(1,1), oma=c(2,2,2,2))
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$maxmonth,zlim=c(1,12),col=turbo(12))
    maps::map(add=TRUE,col='grey',fill=TRUE)
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

pdf(paste0('~/dropbox/working/chlorophyll/plots/gevd_parameter_maps',ocyn,'.pdf'),height=4.5,width=14)
par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location,zlim=c(0,4),col=cols)
    contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    maps::map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2)
    box(lwd=3)
    mtext(expression('a) Location'~mu),adj=0,line=0.5)
    mtext(expression('[mgChl/m'^3*']'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale,zlim=c(-0.5,2),col=cols)
    contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    maps::map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2,labels=NA)
    box(lwd=3)
    mtext(expression('b) Scale'~sigma),adj=0,line=0.5)
    mtext(expression('[mgChl/m'^3*']'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape,zlim=c(0,2.5),col=cols)
    contour(x=lons[lonsi],y=lats[latsi],bathy_atl,add=TRUE,levels=c(50,500,1000,1500,2000,2500,3000))
    maps::map(add=TRUE,col='grey',fill=TRUE)
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

pdf(paste0('~/dropbox/working/chlorophyll/plots/gevd_parameter_se_maps_revision',ocyn,'.pdf'),height=4.5,width=14)
par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))

  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_se/DAT$location,zlim=c(0,0.5),col=cols)
    maps::map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2)
    box(lwd=3)
    mtext(expression('a) sd('*mu*')/'*mu),adj=0,line=0.5)
    mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_se/DAT$scale,zlim=c(0,2),col=cols)
    maps::map(add=TRUE,col='grey',fill=TRUE)
    axis(side=1); axis(side=2,labels=NA)
    box(lwd=3)
    mtext(expression('b) sd('*sigma*')/'*sigma),adj=0,line=0.5)
    mtext(expression('[unitless]'),adj=1.1,line=0.25)
  
  image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_se/DAT$shape,zlim=c(0,0.5),col=cols)
    maps::map(add=TRUE,col='grey',fill=TRUE)
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
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('a) sd('*mu*')/'*mu),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_se/DAT$scale,zlim=c(0,2),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('b) sd('*sigma*')/'*sigma),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_se/DAT$shape,zlim=c(0,0.5),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  axis(side=1); axis(side=2,labels=NA)
  mtext(expression('c) sd('*xi*')/'*xi),adj=0,line=0.5)
  mtext(expression('[unitless]'),adj=1.1,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_half_se/DAT$location_half,zlim=c(0,0.5),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2)
  box(lwd=3)
  mtext(expression('d)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_half_se/DAT$scale_half,zlim=c(0,2),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  #mtext('Halving Sample Size',line=0.5)
  mtext(expression('e)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_half_se/DAT$shape_half,zlim=c(0,0.5),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('f)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_doub_se/DAT$location_doub,zlim=c(0,0.5),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('g)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_doub_se/DAT$scale_doub,zlim=c(0,2),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  #mtext('Doubling Sample Size',line=0.5)
  mtext(expression('h)'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_doub_se/DAT$shape_doub,zlim=c(0,0.5),col=cols)
  maps::map(add=TRUE,col='grey',fill=TRUE)
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
X   <- cbind(c(DAT$qq[,,1,]),c(DAT$qq[,,2,]))
X2  <- X[X[,1]<10 & X[,1]>1E-5 & X[,2]<10 & X[,2]>1E-5,]  #remote extreme outliers where the parameter optimization failed
XX  <- X2[sample(size=5000,1:nrow(X2)),]             #take a random sample for plotting purposes

cor(X2[,1],X2[,2],use='pairwise.complete.obs') #evaluate correlation
cor(log10(X2[,1]),log10(X2[,2]),use='pairwise.complete.obs')

pdf(paste0('~/dropbox/working/chlorophyll/plots/qqplot',ocyn,'.pdf',height=4.5, width=9))
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
pdf(paste0('~/dropbox/working/chlorophyll/plots/qqmap',ocyn,'.pdf',height=5,width=9))
par(mfrow=c(1,1))
  image.plot(x=lons[lonsi],y=lats[latsi],DAT$cormap,col=turbo(20),zlim=c(0.7,1.0),xlab='',ylab='')
  maps::map(add=TRUE,col='grey',fill=TRUE)
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
  maps::map(add=TRUE)
  mtext(expression('a) '*Delta*'BIC Location'~mu~'[mgChl/m'^3*']'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],DAT$chl_trend,zlim=c(-0.1,0.1),col=redblue(20))
  maps::map(add=TRUE)
  mtext(expression('d) Chlorophyll trend [mgChl/m'^3*'/yr]'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],dBIC_scl,zlim=c(0,10),col=turbo(20))
  maps::map(add=TRUE)
  mtext(expression('b) '*Delta*'BIC Scale'~sigma~'[mgChl/m'^3*']'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],DAT$sst_trend,zlim=c(-0.2,0.2),col=redblue(20))
  maps::map(add=TRUE)
  mtext(expression('e) Sea Surface Temperature trend ['*degree*'C/yr]'), adj=0,line=0.25)
  box(lwd=2)
image.plot(x=lons[lonsi],y=lats[latsi],dBIC_shp,zlim=c(0,10),col=turbo(20))
  maps::map(add=TRUE)
  mtext(expression('c) '*Delta*'BIC Shape'~xi~'[unitless]'), adj=0,line=0.25)
  box(lwd=2)
mtext('Latitude',side=2,outer=TRUE,line=0.5)
mtext('Longitude',side=1,outer=TRUE)
dev.off()


##########################################################################
##--MAP OF NONSTATIONARY PARAMETERS--#####################################
##########################################################################
redblue <- colorRampPalette(c('dark blue','blue','white','red','dark red'))

pdf('~/dropbox/working/chlorophyll/plots/gevd_trends_location_scale_shape.pdf',height=4.5,width=14)
par(mfrow=c(1,3),mar=c(1,1,2,4),oma=c(4,4,2,2))

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$location_trnd,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('a) Location'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$scale_trnd,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('b) Scale'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=DAT$shape_trnd,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  axis(side=1); axis(side=2,labels=NA)
  mtext(expression('b) Scale'),adj=0,line=0.5)
  
  mtext('Latitude',outer=TRUE,side=2,line=1.5,cex=1.2)
  mtext('Longitude',outer=TRUE,side=1,line=1.5,cex=1.2)
dev.off()  
  
  
#####################################################################################################
##--MAP OF NONSTATIONARY PARAMETERS WHERE SIGNIFICANT AND WITH CORRELATIONS--########################
#####################################################################################################
loc_sig <- DAT$location_trnd
loc_sig[dBIC_loc<2] <- NA
scl_sig <- DAT$scale_trnd
scl_sig[dBIC_scl<2] <- NA
shp_sig <- DAT$shape_trnd
shp_sig[dBIC_shp<2] <- NA


par(mfcol=c(3,3),mar=c(1,1,2,4),oma=c(4,4,2,2))
image.plot2(x=lons[lonsi],y=lats[latsi],z=loc_sig,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  mtext(expression('a) Location'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=scl_sig,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  mtext(expression('b) Scale'),adj=0,line=0.5)

image.plot2(x=lons[lonsi],y=lats[latsi],z=shp_sig,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  box(lwd=3)
  axis(side=1); axis(side=2,labels=NA)
  mtext(expression('b) Shape'),adj=0,line=0.5)

  mtext('Latitude',outer=TRUE,side=2,line=1.5,cex=1.2)
  mtext('Longitude',outer=TRUE,side=1,line=1.5,cex=1.2)


##################################################################################################
##--SCATTERPLOTS PARAMETER RATE OF CHANGE VS CHL AND SST RATE OF CHANGE--#########################
##################################################################################################
pdf('~/dropbox/working/chlorophyll/plots/scatter_parameter_trends.pdf',height=6,width=9)
par(mfrow=c(2,3),mar=c(2,2,4,2),oma=c(2,2,2,2))
plot(c(DAT$chl_trend),c(loc_sig),pch=19,cex=0.1)#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
  #abline(lm(c(loc_sig) ~ c(DAT$chl_trend)))
  mtext(adj=0.1,line=-2,bquote('r ='~.(round(cor(c(DAT$chl_trend),c(loc_sig),use='pairwise.complete.obs'),3))))  
  mtext(side=2,line=2.5,'Location Parameter Trend')
    
plot(c(DAT$chl_trend),c(scl_sig),ylim=c(-2,2),pch=19,cex=0.1)#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
  mtext(adj=0.1,line=-2,bquote('r ='~.(round(cor(c(DAT$chl_trend),c(scl_sig),use='pairwise.complete.obs'),3))))  
  mtext(side=1,expression('Chlorophyll trend [mgChl/m'^3*'/yr]'),line=3)
  mtext(side=2,line=2.5,'Scale Parameter Trend')
  
plot(c(DAT$chl_trend),c(shp_sig),ylim=c(-2,2),pch=19,cex=0.1)#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
  mtext(adj=0.1,line=-2,bquote('r ='~.(round(cor(c(DAT$chl_trend),c(shp_sig),use='pairwise.complete.obs'),3))))  
  mtext(side=2,line=2.5,'Shape Parameter Trend')
  
plot(c(DAT$sst_trend),c(loc_sig),pch=19,cex=0.1)#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
  mtext(adj=0.1,line=-2,bquote('r ='~.(round(cor(c(DAT$sst_trend),c(loc_sig),use='pairwise.complete.obs'),3))))  
  mtext(side=2,line=2.5,'Location Parameter Trend')
  
plot(c(DAT$sst_trend),c(scl_sig),ylim=c(-2,2),pch=19,cex=0.1)#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
  mtext(adj=0.1,line=-2,bquote('r ='~.(round(cor(c(DAT$sst_trend),c(scl_sig),use='pairwise.complete.obs'),3))))  
  mtext(side=2,line=2.5,'Scale Parameter Trend')
  mtext(side=1,expression('Sea Surface Temperature trend ['*degree*'C/yr]'),line=3)
  
plot(c(DAT$sst_trend),c(shp_sig),ylim=c(-2,2),pch=19,cex=0.1)#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
  mtext(adj=0.1,line=-2,bquote('r ='~.(round(cor(c(DAT$sst_trend),c(shp_sig),use='pairwise.complete.obs'),3))))  
  mtext(side=2,line=2.5,'Shape Parameter Trend')
dev.off()  
  

##--CALCULATE CORRELATIONS BETWEEN ENVIRONMENTAL AND PARAMETER TRENDS--########### 
cor(c(DAT$chl_trend),c(loc_sig),use='pairwise.complete.obs')#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
cor(c(DAT$chl_trend),c(scl_sig),use='pairwise.complete.obs')#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
cor(c(DAT$chl_trend),c(shp_sig),use='pairwise.complete.obs')#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  

cor(c(DAT$sst_trend),c(loc_sig),use='pairwise.complete.obs')#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
cor(c(DAT$sst_trend),c(scl_sig),use='pairwise.complete.obs')#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  
cor(c(DAT$sst_trend),c(shp_sig),use='pairwise.complete.obs')#,xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))  




################################################################################################
##--MAPS PARAMETER RATES OF CHANGE VS RATES OF CHANGE IN CHL AND SST--##########################
################################################################################################
pdf('~/dropbox/working/chlorophyll/plots/nonstationary_delta_theta_delta_chl_sst.pdf',height=7,width=8)
par(mfrow=c(3,2),mar=c(1,1,2,4),oma=c(4,4,2,2))
  image.plot2(x=lons[lonsi],y=lats[latsi],z=loc_sig,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2)
  box(lwd=3)
  mtext(expression('a) Location trend [mgChl/m'^3*'/yr]'),adj=0,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],DAT$chl_trend,zlim=c(-0.1,0.1),col=redblue(20))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2,labels=NA)
  mtext(expression('d) Chlorophyll trend [mgChl/m'^3*'/yr]'), adj=0,line=0.25)
  box(lwd=2)

image.plot2(x=lons[lonsi],y=lats[latsi],z=scl_sig,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1,labels=NA); axis(side=2)
  box(lwd=3)
  mtext(expression('b) Scale trend [mgChl/m'^3*'/yr]'),adj=0,line=0.25)

image.plot2(x=lons[lonsi],y=lats[latsi],DAT$sst_trend,zlim=c(-0.2,0.2),col=redblue(20))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2,labels=NA)
  mtext(expression('e) Sea Surface Temperature trend ['*degree*'C/yr]'), adj=0,line=0.25)
  box(lwd=2)

image.plot2(x=lons[lonsi],y=lats[latsi],z=shp_sig,zlim=c(-0.25,0.25),col=redblue(12))
  maps::map(add=TRUE,col='grey',fill=TRUE)
  axis(side=1); axis(side=2)
  box(lwd=3)
  axis(side=1); axis(side=2,labels=NA)
  mtext(expression('c) Shape trend [/yr]'),adj=0,line=0.25)

  mtext(side=1,'Longitude',line=1.0,outer=TRUE)
  mtext(side=2,'Latitude',line=1.5,outer=TRUE)
dev.off()




