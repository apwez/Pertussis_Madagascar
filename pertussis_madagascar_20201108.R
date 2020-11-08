library(tidyverse)
library(lubridate)
library(beeswarm)
library(RColorBrewer)


#### 2. Focal data: pertussis serology  #########
df1 <- read.csv("pertussis_sero_data.csv")

df1$titer_group <- cut(df1$Titer..IU.ml., breaks = c(-1,5,40,100, 10000), labels = c('<5', '5-40', '40-100', '>=100'))
df1$age_group <- cut((df1$Age..years.*12 + df1$Age..months.)/12, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), labels = c('6m-11m', '1y', '2y', '3y', '4y', '5y', '6y', '7y', '8y', '9y', '10y', '11y', '12y', '13y', '14-15y'))

### Clean up variables
df1$Site.d.étude <- trimws(df1$Site.d.étude)
df1$Immunization.card <- trimws(df1$Immunization.card)
df1$Date.Sampling <- as.Date(df1$Date.Sampling, "%d/%m/%Y")
df1$Date.of.Birth <- as.Date(df1$Date.of.Birth, "%d/%m/%Y")
df1$VPO <- as.Date(df1$VPO, "%d/%m/%Y")
df1$VP1 <- as.Date(df1$VP1, "%d/%m/%Y")
df1$VP2 <- as.Date(df1$VP2, "%d/%m/%Y")
df1$VP3 <- as.Date(df1$VP3, "%d/%m/%Y")
#df1$VP4 <- as.Date(df1$VP4, "%d/%m/%Y")
vaccinated <- 1*(!is.na(df1$VPO) & !is.na(df1$VP1) & !is.na(df1$VP2) & !is.na(df1$VP3))
vaccinated3 <- 1*(!is.na(df1$VPO) & !is.na(df1$VP1) & !is.na(df1$VP2))
vaccinated2 <- 1*(!is.na(df1$VPO) & !is.na(df1$VP1))
vaccinated1 <- 1*(!is.na(df1$VPO))

### age 
df1$Age <- as.numeric(df1$Date.Sampling-df1$Date.of.Birth)/365
df1$Age[is.na(df1$Age)] <- df1$Age..years.[is.na(df1$Age)]+sample(-182:182,size=1,replace=TRUE)/365
df1$Age[df1$Age<1] <- (df1$Age..months.[df1$Age<1]-1)/12

### Clean up district variables
sites.numeric <- as.numeric(as.factor(df1$District.N.))
names.sites <- rep(NA,max(sites.numeric,na.rm=TRUE))
for (j in 1:length(names.sites)) names.sites[j] <- df1$District.N.[sites.numeric==j][1]
names.sites.long <- c("Antananarivo", "Midongy Atsimo", "Mahajanga I", "Antsalova", "Toliara I")


### run code in STAN 
require(lubridate)
source('~/Dropbox (Princeton)/Madagascar-Rubella-IgG/source/additionalExperiments.R', chdir = TRUE)  #where making data for stan resides
date.birth <- decimal_date(df1$Date.of.Birth)
date.collection <- decimal_date(df1$Date.Sampling)
date.birth[is.na(date.birth)] <- date.collection[is.na(date.birth)]-df1$Age[is.na(date.birth)]  ##note now using randomization from above

mat <- cbind(date.birth-2001,date.collection-2001,df1$Age, sites.numeric)
seroneg.ind.here <- mat[which(df1$Titer..IU.ml.<2,arr.ind=TRUE),]
seropos.ind.here <-  mat[which(df1$Titer..IU.ml.>=2,arr.ind=TRUE),]


## 
##' Extract for stan 
##'
##' @param seroneg.ind  - matrix with one row per person, year birth, year measure, age 
##'                       these are adjusted to make the origin year 0, so some years birth are negative       
##' @param seropos.ind  - same but for seropositive  
##' @param duration.maternal.immunity - how long is maternal immunity expected to last  
##' @param start.year - when in the year is trough 
##' @param min.year - which is earliest year to start with (only want to do last 15 for ex)   
##' @param imposeT - if you want the number of columns set (to make cases match the seropos)
##' ...
##'
##' @return list containing x[N,T] fraction each year individual present
##'                         y[N];    //if seropositive or seronegative for each ind

extractForStan <- function(seroneg.ind,seropos.ind, 
                           duration.maternal.immunity=9/12, 
                           start.year=0.5,min.year=35, imposeT=NULL) {
  
  ## restrict to focal years (otherwise trying too widely)
  ## negative values tell how many years alive must account for before start focal
  seropos.ind <- seropos.ind-min.year+start.year  
  seroneg.ind <- seroneg.ind-min.year+start.year   ## !add start.year here to ease line up with tot.year
  
  ## pull out any individuals who lived BEFORE the focal span
  if (!is.null(nrow(seropos.ind))) seropos.ind <- seropos.ind[seropos.ind[,2]>0, ]
  if (!is.null(nrow(seroneg.ind))) seroneg.ind <- seroneg.ind[seroneg.ind[,2]>0, ]
  
  ## adjust for maternal immunity 
  if (!is.null(nrow(seroneg.ind))) { 
    seroneg.ind[,1] <- seroneg.ind[,1]+duration.maternal.immunity
  } else {
    seroneg.ind[1] <- seroneg.ind[1]+duration.maternal.immunity}
  if (!is.null(nrow(seropos.ind))) {
    seropos.ind[,1] <- seropos.ind[,1]+duration.maternal.immunity
  } else {
    seropos.ind[1] <- seropos.ind[1]+duration.maternal.immunity}
  
  ## those that contain no values,  become NAs here - replace with numeric 0
  if (is.na(sum(seroneg.ind))) seroneg.ind <- c()
  if (is.na(sum(seropos.ind))) seropos.ind <- c()
  
  ##make the response and the matrix
  if (!is.null(nrow(seropos.ind)) &!is.null(nrow(seroneg.ind))) {
    yVec <- c(rep(0,nrow(seroneg.ind)),rep(1,nrow(seropos.ind)))
    if (!is.null(imposeT)) T <- imposeT else T <- floor(max(c(seroneg.ind[,2],seropos.ind[,2])))+1
    xMat <- matrix(0,nrow(seroneg.ind)+nrow(seropos.ind),T)}
  
  if (!is.null(nrow(seropos.ind)) & is.null(nrow(seroneg.ind))) {
    yVec <- c(rep(0,1),rep(1,nrow(seropos.ind)))
    if (!is.null(imposeT)) T <- imposeT else T <- floor(max(c(1,seropos.ind[,2])))+1
    xMat <- matrix(0,1+nrow(seropos.ind),T)}
  
  if (is.null(nrow(seropos.ind)) & !is.null(nrow(seroneg.ind))){ 
    yVec <- c(rep(0,nrow(seroneg.ind)),rep(1,1))
    if (!is.null(imposeT)) T <- imposeT else T <- floor(max(c(seroneg.ind[,2],1)))+1
    xMat <- matrix(0,nrow(seroneg.ind)+1,T)}
  
  if (is.null(nrow(seropos.ind)) & is.null(nrow(seroneg.ind))){ 
    yVec <- c(rep(0,1),rep(1,1))
    if (!is.null(imposeT)) T <- imposeT else T <- floor(max(c(seroneg.ind,1)))+1
    xMat <- matrix(0,2,T)}
  
  #print(dim(xMat))
  #print(seroneg.ind)
  
  if (!is.null(nrow(seroneg.ind))) {
    for (kk in 1:nrow(seroneg.ind)) {  ##loop over individuals
      #put 1s in all years that exist in, with index starting in 2
      xMat[kk,(max(floor(seroneg.ind[kk,1]),0)+1):(floor(seroneg.ind[kk,2])+1)] <- 1
      
      #in first year, put number of years alive before sample begins
      xMat[kk,1] <- (seroneg.ind[kk,1]<0)*abs(seroneg.ind[kk,1])
      
      ## modulate margins by how much of the year is included
      # in last year, reduce it by fraction of year that present
      xMat[kk,(floor(seroneg.ind[kk,2])+1)] <- xMat[kk,(floor(seroneg.ind[kk,2])+1)]*(seroneg.ind[kk,2]-floor(seroneg.ind[kk,2]))
      
      # if first year>0, reduce it by time of birth in year and maternal immunity
      ## xxj ## if (seroneg.ind[kk,1]>=0) xMat[kk,(floor(seroneg.ind[kk,1])+1)] <- xMat[kk,(floor(seroneg.ind[kk,1])+1)]*(1-(seroneg.ind[kk,1]-floor(seroneg.ind[kk,1])))
      if (seroneg.ind[kk,1]>=0) xMat[kk,(floor(seroneg.ind[kk,1])+1)] <- (1-(seroneg.ind[kk,1]-floor(seroneg.ind[kk,1])))
      
      n.neg <- nrow(seroneg.ind)
      
    }} else {
      if(!is.null(seroneg.ind)){
        xMat[1,(max(floor(seroneg.ind[1]),0)+1):(floor(seroneg.ind[2])+1)] <- 1
        xMat[1,1] <- (seroneg.ind[1]<0)*abs(seroneg.ind[1])
        xMat[1,(floor(seroneg.ind[2])+1)] <- xMat[1,(floor(seroneg.ind[2])+1)]*(seroneg.ind[2]-floor(seroneg.ind[2]))
        if (seroneg.ind[1]>=0) xMat[1,(floor(seroneg.ind[1])+1)] <- xMat[1,(floor(seroneg.ind[1])+1)]*(1-(seroneg.ind[1]-floor(seroneg.ind[1])))
        n.neg <- 1 } else {n.neg <- 0}
    }
  # print(dim(xMat))
  if (!is.null(nrow(seropos.ind))) {
    
    for (kk in 1:nrow(seropos.ind)) {  ##loop over individuals
      #print(n.neg+kk)
      #print((max(floor(seropos.ind[kk,1]),0)+1):(floor(seropos.ind[kk,2])+1))
      #put 1s in all years that exist in, with index starting in 2
      xMat[n.neg+kk,(max(floor(seropos.ind[kk,1]),0)+1):(floor(seropos.ind[kk,2])+1)] <- 1
      #print("here")
      #in first year, put number of years alive before sample begins
      xMat[n.neg+kk,1] <- (seropos.ind[kk,1]<0)*abs(seropos.ind[kk,1])
      
      ## modulate margins by how much of the year is included
      # in last year, reduce it by fraction of year that present
      xMat[n.neg+kk,(floor(seropos.ind[kk,2])+1)] <- xMat[n.neg+kk,(floor(seropos.ind[kk,2])+1)]*(seropos.ind[kk,2]-floor(seropos.ind[kk,2]))
      # if first year>0, reduce it by time of birth in year and maternal immunity
      if (seropos.ind[kk,1]>=0) xMat[n.neg+kk,(floor(seropos.ind[kk,1])+1)] <- xMat[n.neg+kk,(floor(seropos.ind[kk,1])+1)]*(1-(seropos.ind[kk,1]-floor(seropos.ind[kk,1])))
      
    }}else {
      if(!is.null(seropos.ind)){
        xMat[n.neg+1,(max(floor(seropos.ind[1]),0)+1):(floor(seropos.ind[2])+1)] <- 1
        xMat[n.neg+1,1] <- (seropos.ind[1]<0)*abs(seropos.ind[1])
        xMat[n.neg+1,(floor(seropos.ind[2])+1)] <- xMat[1,(floor(seropos.ind[2])+1)]*(seropos.ind[2]-floor(seropos.ind[2]))
        if (seropos.ind[1]>=0) xMat[1,(floor(seropos.ind[1])+1)] <- xMat[1,(floor(seropos.ind[1])+1)]*(1-(seropos.ind[1]-floor(seropos.ind[1])))
      }}
  
  return(list(xMat=xMat,yVec=yVec))
  
}



a1 <- extractForStan(seroneg.ind=seroneg.ind.here[,1:3],seropos.ind=seropos.ind.here[,1:3],
                     duration.maternal.immunity=1.5/12, 
                     start.year=0.5,min.year=1, imposeT=NULL)
a1$ages <- c(seroneg.ind.here[,3],seropos.ind.here[,3])
a1$locs <- c(seroneg.ind.here[,4],seropos.ind.here[,4])


### run stan code to fit FOI
require(rstan)

fit1d <- stan(file="stan1d.stan", 
              data=list(x=a1$xMat[,1:15],y=a1$yVec,N=nrow(a1$xMat),T=ncol(a1$xMat)-1, 
              lambdaStart=0.0001,lambdaStartSigma=1, ## note currently not used
              vStart=rep(0.5,length(unique(a1$locs))),vStartSigma=rep(5,length(unique(a1$locs))), 
              nloc=length(unique(a1$locs)),
              locs=a1$locs), chains=4, iter=4000)
rc1d <- rstan::extract(fit1d)
exp(median(extract_log_lik(fit1d,parameter_name = "logLik")))








