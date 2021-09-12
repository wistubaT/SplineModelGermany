library('readxl')

effectiveIFR <- function(ifrSrc){
  # Read data of cases per age group and calendar week in Germany and prepare the data frame
  dat <- read_excel('data/cases_per_age-group_calendar-week.xlsx', sheet = 2)
  dat <- dat[-1,]
  dat <- data.frame(dat)
  row.names(dat) <- dat[,1]
  dat <- dat[seq(dim(dat)[1],1),]
  dat <- dat[,-1]
  
  # Read data of population per age in Germany and prepare the data frame
  pop <- read.csv('data/population.csv', sep = ';', row.names = 1)
  pop[17,] <- colSums(pop)
  row.names(pop)[17] <- 'Deutschland'
  pop <- pop*1000
  pop_ant <- pop
  
  # Extract the integers of the different age groups in dat
  ags <- matrix(rep(0,2*dim(dat)[1]), ncol = 2)
  for(i in 1:(dim(dat)[1]-1)){
    ag <- strsplit(rownames(dat)[i], '-')
    ags[i,1] <- as.integer(ag[[1]][1])+1
    ags[i,2] <- as.integer(ag[[1]][2])+1
  }
  ags[i+1,1] <- ags[i,2]+1
  ags[i+1,2] <- dim(pop)[2]
  
  # Read the data with the desired ifr estimates
  ifr <- read.csv(paste0('data/ifr_', tolower(ifrSrc),'.csv'), sep = ';')
  ifr[1,1] <- 0
  ifr[dim(ifr)[1],2] <- 100
  ifr[,c(1,2)] <- ifr[,c(1,2)]+1
  ifrvec <- NULL
  for(i in 1:dim(ifr)[1]){
    ifrvec <- c(ifrvec, rep(ifr[i,3], ifr[i,2]-ifr[i,1]+1))
  }
  
  # Initialize matrices for cases per age group and calendar week and ifr per calendar week
  cpagcw <- matrix(rep(0, dim(dat)[1]*dim(dat)[2]), nrow = dim(dat)[1])
  ifrpcw <- matrix(rep(0,dim(pop)[1]*dim(dat)[2]), nrow = dim(pop)[1])
  
  # Estimate the cases per age group and calendar week in Germany
  
  for(i in 1:dim(dat)[1]){
    popag <- sum(pop[17,ags[i,1]:ags[i,2]])
    for(j in 1:dim(dat)[2]){
      cpagcw[i,j] <- dat[i,j] / popag
    }
  }
  
  # Estimate the ifr per calendar week as a weighted harmonic mean
  
  for(k in 1:dim(pop)[1]){
    for(j in 1:dim(cpagcw)[2]){
      infest <- NULL
      for(i in 1:dim(cpagcw)[1]){
        infest <- c(infest, cpagcw[i,j]*as.matrix(pop[k,ags[i,1]:ags[i,2]]))
      }
      ifrpcw[k,j] <- (infest %*% ifrvec) / sum(infest)
    }
  }
  return(ifrpcw)
}