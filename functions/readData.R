read_covid19_data_rki <- function(){
  # Download new case data if it is available
  if(!file.exists('last_download_date.Rdata')){
    last_download_date <- today() - 1
    save(last_download_date, file = 'last_download_date.Rdata')
  }
  load('last_download_date.Rdata')
  if(last_download_date != today()){
    print('New case data is available. Downloading new data...')
    download.file('https://www.arcgis.com/sharing/rest/content/items/f10774f1c63e40168479a1feb6c7ca74/data','data/RKI_COVID19.csv')
    last_download_date <- today()
    save(last_download_date, file = 'last_download_date.Rdata')
  }
  
  # Read case data
  dat <- read.csv('data/RKI_COVID19.csv', sep = ',', stringsAsFactors = FALSE)
  dat <- data.frame(dat, stringsAsFactors = FALSE)
  dat <- dat[which(dat$NeuerFall %in% c(0,1)),]
  dat <- data.frame(numCases = dat$AnzahlFall,
                  federalState = dat$Bundesland,
                  repDate = dat$Meldedatum, stringsAsFactors = FALSE)
  dat$repDate <- strtrim(dat$repDate,10)
  dat$repDate <- as.Date(dat$repDate, format = '%Y/%m/%d')
  dat$federalState[which(dat$federalState=='Baden-Württemberg')] <- 'Baden-Wuerttemberg'
  dat$federalState[which(dat$federalState=='Thüringen')] <- 'Thueringen'
  
  # Read death data
  deaths <- read.csv('data/RKI_deaths.csv', sep = ';', row.names = 1, stringsAsFactors = FALSE)
  colnames(deaths) <- c('Baden-Wuerttemberg','Bayern','Berlin','Brandenburg','Bremen','Hamburg','Hessen',
                      'Mecklenburg-Vorpommern','Niedersachsen','Nordrhein-Westfalen','Rheinland-Pfalz',
                      'Saarland','Sachsen','Sachsen-Anhalt','Schleswig-Holstein','Thueringen')
  row.names(deaths) <- as.Date(row.names(deaths), format = '%d.%m.%Y')
  
  if(max(as.Date(row.names(deaths))) != max(dat$repDate)){
    print('Both datasets have differing max dates. Adjusting both datasets...')
    maxDate <- min(max(as.Date(row.names(deaths))), max(dat$repDate))
    print(paste('The last day with observed data was adjusted to be: ', maxDate))
    deaths <- deaths[which(as.Date(row.names(deaths)) <= maxDate),]
    dat <- dat[which(dat$repDate<=maxDate),]
  }

  regions <- data.frame(Regions = unique(dat$federalState)[order(unique(dat$federalState))], stringsAsFactors = FALSE)
  n <- length(seq.Date(from = as.Date('2019-12-30'), to = max(dat$repDate), by = 1))
  d <- data.frame(DateRep = rep(seq.Date(from = as.Date('2019-12-30'), to = max(dat$repDate), by = 1), 16),
                Cases = rep(0,n*16),
                Deaths = rep(0,n*16),
                Region = rep(regions$Regions,n)[order(rep(regions$Regions,n))],
                stringsAsFactors = FALSE)
  
  # Count cases per day and federal state
  for(i in 1:dim(d)[1]){
    cases_per_day_fedState <- dat$numCases[intersect(which(dat$repDate==d$DateRep[i]),which(dat$federalState==d$Region[i]))]
    d$Cases[i] <- sum(cases_per_day_fedState)
  }
  
  cases <- data.frame(rep(0,n))
  k <- 1
  for(reg in regions$Regions){
    cases[,k] <- d$Cases[which(d$Region==reg)]
    k <- k+1
  }
  colnames(cases) <- regions$Regions
  rownames(cases) <- d$DateRep[which(d$Region==reg)]
  d$Deaths <- unlist(deaths)
  return(list(d = d, cases = cases, deaths = deaths))
}
