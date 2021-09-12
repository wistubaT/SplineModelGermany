addGermany <- function(d, cases, deaths){
  casesD <- data.frame(Cases = rowSums(cases))
  deathsD <- data.frame(Deaths = rowSums(deaths))
  dates <- data.frame(DateRep = unique(d$DateRep)[order(unique(d$DateRep))])
  dD <- cbind(dates, casesD, deathsD, Region = 'Deutschland')
  rownames(dD) <- seq(1, dim(dates)[1])
  
  d2 <- rbind(d, dD)
  cases2 <- cbind(cases,casesD)
  colnames(cases2)[17] <- 'Deutschland'
  deaths2 <- cbind(deaths,deathsD)
  colnames(deaths2)[17] <- 'Deutschland'
  
  return(list(d = d2, cases = cases2, deaths = deaths2))
}
