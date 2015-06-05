
#little function to calculate VPD
getVPD<-function(Ta,RH){  
  svp <- 0.6108 * exp(17.27 * Ta / (Ta + 237.3)) #caclulate saturation vapor pressure, in kPa
  avp <- RH*svp/100
  VPD = svp-avp
  return(VPD)
}

# Function for converting all columns in a dataframe to numeric
numericdfr <- function(dfr){
  for(i in 1:ncol(dfr)){
    options(warn=-1)
    oldna <- sum(is.na(dfr[,i]))
    num <- as.numeric(as.character(dfr[,i]))
    if(sum(is.na(num)) > oldna)
      next
    else
      dfr[,i] <- num
  }
  options(warn=0)
  return(dfr)
}

standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#----------------------------------------------------------------------
# Round POSIXct date/time down to a number of time-levels.
#
round.POSIXct <- function(x, units = c("mins", "5 mins", "10 mins", "15 mins", "quarter hours", "30 mins", "half hours", "hours","6-hours")){
  #if(is.numeric(units)) units <- as.character(units)
  #units <- match.arg(units)
  r <- switch(units,
              "mins" = 60,
              "5 mins" = 60*5,
              "10 mins" = 60*10,
              "15 mins"=, "quarter hours" = 60*15,
              "30 mins"=, "half hours" = 60*30,
              "hours" = 60*60,
              "4 hours" = 60*60*4,
              "6 hours" = 60*60*6,
              "day" = 60*60*24)
  H <- as.integer(format(x, "%H"))
  M <- as.integer(format(x, "%M"))
  S <- as.integer(format(x, "%S"))
  D <- format(x, "%Y-%m-%d",tz = "GMT")
  secs <- 3600*H + 60*M + S
  as.POSIXct(round(secs/r)*r, origin=D,tz="GMT")
}


#return saturating vapour pressure based on air temperature, in millibars
get.es <- function(temp){
  es <- 6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + temp)))
  return(es)
}
