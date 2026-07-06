# The following functions are taken/adapted from the code repository accompanying:
# Tebaldi et al. (2021) Nature Climate Change
# https://github.com/DOE-ICoM/tebaldi-etal_2021_natclimchange

# Modifications: Khin Nawarat (2026, @IHE Delft)
# Modifications in rasmussenScenarios (removing 2 GWL scenarios (2°C+ and 5°C+) representing alternative SEJ-based high-end ice-sheet melt assumptions)

library(matlib)
library(MASS)
library(extRemes)
library(class)
library(R.matlab)
library(RColorBrewer)
library(maps)

compute.variance.covariance<-function(scale,shape,lambda,m){
    temp<-matrix(c(2*scale^2,scale,scale,1+shape),2,2)
    N<-lambda*m
    frac<-(1+shape)/N
    return(frac*temp)}



simulate.bivariate<-function(n,scale,shape,lambda,m){
    mu<-c(scale,shape)
    Sigma<-compute.variance.covariance(scale,shape,lambda,m)
    if(det(Sigma)>0)
    return(mvrnorm(n,mu,Sigma))
    else return(matrix(NA,n,length(mu)))
}
RL.N<-function(N,shape,scale,thresh,lambda){thresh+(scale/shape)*((N*lambda)^shape-1)}   #derive return levels based on GP parameters
RP<-function(Z,shape,scale,thresh,lambda){
    rps<-ifelse(Z<=thresh,1,(1/(lambda))*(((shape/scale)*(Z-thresh)+1)^(1/shape)))
    rps[is.na(rps)]<-1000
    rps[rps<1]<-1
    rps[rps>1000]<-1000
    rps}

find.thresh<-function(Z,shape,scale,lambda,RP=100){
    thresh=Z-(scale/shape)*((lambda*RP)^shape-1)
    return(thresh)
}

unfold.quantiles<-function(qvalues,qs=c(1,5,16.7,50,83.3,95,99,99.5,99.9),n=1000){
  
  if(any(qs>1))qs<-qs/100
  samplesizes<-diff(qs)*n
  
  qsample<-numeric(0)
  for(ss in seq(length(samplesizes))){
    a<-qvalues[ss]
    b<-qvalues[ss+1]
    qsample<-c(qsample,runif(samplesizes[ss],a,b))
  }
  return(qsample)
}



g.coord<-function(x,y){
    pp<-par("usr")
    xr<-pp[2]-pp[1]
    yr<-pp[4]-pp[3]
    xpos<-x*xr+pp[1]
    ypos<-y*yr+pp[3]
    return(list(x=xpos,y=ypos))
}


"%&%"<-function(x,y)paste0(x,y)

find.match<-function(match1,match2,k=1,thresh=1){ #match2 is the larger set
    indy<-as.numeric(knn(train=cbind(match2$longitude,match2$latitude),test=cbind(match1$longitude,match1$latitude),cl=seq(nrow(match2)),k=k))#the order to apply to match1 to match match2
    distances<-numeric(0)
    for(i in 1:nrow(match1))distances[i]<-sqrt((match1$longitude[i]-match2$longitude[indy][i])^2+(match1$latitude[i]-match2$latitude[indy][i])^2)
    indyclean<-indy[distances<=1]
    rowclean<-seq(nrow(match1))[distances<=1]
    temp<-cbind(rowclean,indyclean)
    dimnames(temp)<-list(NULL,c("match1","match2"))
    return(temp)}

taper<-function(x,q1=0.05,q2=0.95){
    qx<-quantile(x,prob=c(q1,q2),na.rm=TRUE)
    x[x<qx[1]]<-qx[1]
    x[x>qx[2]]<-qx[2]
    return(x)}


get.early<-function(x,y){
    ux<-unique(x)
    uy<-numeric(length(ux))
    for(i in seq(length(uy))){
        uy[i]<-min(y[x==ux[i]])
    }
    return(uy)
}

transparent_col <- function(color, percent = 50, name = NULL) {
    #      color = color name
    #    percent = % transparency
    #       name = an optional name for the color

    ## Get RGB values for named color
    rgb.val <- col2rgb(color)

    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)

    ## Save the color
    invisible(t.col)
}


rasmussenScenarios<-c("1p5degree","2p0degree","2p5degree","3p0degree","4p0degree","5p0degree")



