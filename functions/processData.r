source('functions/effectiveIFR.R')

processData <- function(region, d, pop,N2, indRegion, ifrSrc){
  ifr.by.bl <- effectiveIFR(ifrSrc)
  serial.interval = read.csv('data/serial_interval.csv')
  
  # Pads serial interval with 0 if N2 is greater than the length of the serial
  # interval array
  if (N2 > length(serial.interval$fit)) {
    pad_serial.interval <- data.frame(
      'X'=(length(serial.interval$fit)+1):N2,
      'fit'=rep(1e-17, max(N2-length(serial.interval$fit), 0 ))
    )
    serial.interval <- rbind(serial.interval, pad_serial.interval)
  }
  
  # various distributions required for modeling
  x1 <- rlnorm(1e6,1.62,0.42) # incubation period from https://github.com/HopkinsIDD/ncov_incubation
  mean <- 17.8; cv <- 0.45 # onset to death
  x2 <- rgammaAlt(1e6,mean,cv) # onset-to-death distribution
  ecdf.saved <- ecdf(x1+x2)
  
  # Convert the weekly ifr estimates to daily ifr estimates (constant for each day in a week)
  IFRtemp <- as.numeric(ifr.by.bl[indRegion,])
  IFR.start <- as.Date('2020-03-02') - 10
  IFR.end <- as.Date('2020-03-02') - 10 + length(IFRtemp)*7 - 1
    
  IFR <- NULL
  for(i in 1:length(IFRtemp)) IFR <- c(IFR, rep(IFRtemp[i], 7))
  IFR <- c(rep(IFR[1], as.numeric(IFR.start-dmy('31/12/2019'))), IFR)
    
  if (N2 > length(IFR)) {
    pad_IFR <- rep(IFR[length(IFR)], length((length(IFR)+1):N2))
    IFR <- c(IFR, pad_IFR)
  }
    
  d <- d[order(as.Date(d$DateRep, format = '%Y-%m-%d')),]  # ensure date ordering
  
  # Determine the beginning of modelling
  index <- which(d$Cases>0)[1]
  index1 <- which(cumsum(d$Deaths)>=10)[1] # also 5
  index2 <- max(c(index1-60,1))
    
  print(sprintf('First non-zero cases is on day %d, and 60 days before 10 deaths is day %d',index,index2))
  d <- d[index2:nrow(d),]
  
  # Pad the ifr data frame according to the adjusted time interval
  IFR <- IFR[index2:length(IFR)]
  pad_IFR <- rep(IFR[length(IFR)], length((length(IFR)+1):N2))
  IFR <- c(IFR, pad_IFR)
    
  dates <- d$DateRep
  N <- length(d$Cases)
  print(sprintf('%s has %d days of data',region,N))
    
  daysDiff <- N2 - N
  if(daysDiff < 0) {
    print(sprintf('%s: %d', Country, N))
    print('ERROR!!!! increasing N2')
    N2 <- N
    daysDiff <- N2 - N
  }
    
  convolution <- function(u) (ecdf.saved(u))
    
  discPi <- rep(0,N2) # discretized distribution function PI
  discPi[1] <- (convolution(1.5) - convolution(0))
  for(i in 2:N2) {
    discPi[i] <- (convolution(i+.5) - convolution(i-.5)) 
  }
  deaths <- c(as.vector(as.numeric(d$Deaths)),rep(-1,daysDiff))
  cases <- c(as.vector(as.numeric(d$Cases)),rep(-1,daysDiff))
  stan_data <- list(N=N, deaths=deaths, discPi=discPi, N0=6, cases=cases, SI=serial.interval$fit[1:N2], 
                    features=NULL, EpidemicStart=index1+1-index2, pop=pop, N2=N2, x=1:N2, 
                    B2=nrow(t(bs(seq(1,N2,1) , knots=seq(1,N2,14), degree=3, intercept = FALSE))), B=t(bs(seq(1,N2,1) , knots=seq(1,N2,14), degree=3, intercept = FALSE)),
                    IFR=c(IFR))
  
  
  return(list('stan_data' = stan_data, 'dates' = dates, 'reported_cases'=as.vector(as.numeric(d$Cases)), 
              'deaths_by_region' = as.vector(as.numeric(d$Deaths))))
}
