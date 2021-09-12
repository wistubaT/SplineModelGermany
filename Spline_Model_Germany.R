setwd(INSERT_FILE_LOCATION)
Sys.setlocale('LC_TIME', 'English') 
library(rstan)
library(lubridate)
library(dplyr)
library(tidyr)
library(EnvStats)
library(gridExtra)
library(splines)
library(ggplot2)
library(bayesplot)
source('functions/processData.r')
source('functions/readData.r')
source('functions/addGermany.R')
source('functions/plotFunctionsFS.R')

#### Initialize model parameters ####

modeltype <- 'HALF'    # Specify HMC-Sampler: DEBUG   (testing for errors in coding)
                       #                      HALF    (testing the model with a low chain and iteration counts)
                       #                      NORMAL  (running the full model)
ifrSrc <- 'BR'         # Select the source of age-specific ifr estimates to be used (BR: Brazeau et al., OD: O'Driscoll at al., LE: Levin et al.)
indRegion <-  2        # Select your federal state of interest (number corresponds to index
                       # alphabetically ordered German Federal states); remark: 17=Germany

#### Prepare Data ####

regionsList <- c('Baden-Wuerttemberg', 'Bayern', 'Berlin', 'Brandenburg', 'Bremen', 'Hamburg', 'Hessen',
                 'Mecklenburg-Vorpommern', 'Niedersachsen', 'Nordrhein-Westfalen', 'Rheinland-Pfalz', 
                 'Saarland', 'Sachsen', 'Sachsen-Anhalt', 'Schleswig-Holstein', 'Thueringen', 'Deutschland')

# Read the case and death data per day and federal state
readdata <- read_covid19_data_rki()
d <- readdata$d
deaths <- readdata$deaths
cases <- readdata$cases

if(indRegion==17){
  source('functions/plotFunctions.R')
  newD <- addGermany(d, cases, deaths)
  d <- newD$d
  deaths <- newD$deaths
  cases <- newD$cases
}
region <- regionsList[indRegion]

d <- d[which(d$Region == region),]
cases <- cases[,which(colnames(cases) == region)]
deaths <- deaths[,which(colnames(deaths) == region)]

N2 <- (max(d$DateRep) - min(d$DateRep) + 1)[[1]] # number of days with available data

# Read the population data
pops <- read.csv('data/population.csv', stringsAsFactors = FALSE, sep = ';')
rownames(pops) <- pops[,1]
pops <- pops[,-1]*1000
pops[17,] <- colSums(pops)
rownames(pops)[17] <- 'Deutschland'
pops <- rowSums(pops)

# Process the data so it can be used in the .stan file
processedData <- process_covariates_each(region = region, d = d , pop = pops[indRegion], 
                                          N2 = N2, indRegion = indRegion, ifrSrc = ifrSrc)
stan_data = processed_data$stan_data
dates = processed_data$dates
deaths_by_country = processed_data$deaths_by_country
reported_cases = processed_data$reported_cases

#### Model Fitting ####

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(file = 'splineModel.stan') # Auswahl des Modells

# Different sampling methods
if(toupper(modeltype)=='DEBUG'){
  fit = sampling(m,data=stan_data,iter=40,warmup=20,chains=2)
}else if(toupper(modeltype)=='HALF'){
  fit = sampling(m,data=stan_data,iter=1000,warmup=500,chains=4,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 15))
}else if(toupper(modeltype)=='NORMAL'){
  fit = sampling(m,data=stan_data,iter=2000,warmup=1000,chains=8,thin=1,control = list(adapt_delta = 0.99, max_treedepth = 15))
}

# Extract the A-posteriori samples
out <- rstan::extract(fit)
setwd('Rdata')
save(readdata, processed_data, m, fit, out, file = paste0('results',paste0(today()), '_', toupper(modeltype), '_', region, '_', floor(runif(1,min=1,max=1e5)),'_', ifrSrc,'.Rdata'))
setwd('..')

#### Prepare Plotting ####
N <- stan_data$N # Number of modelled days

predicted_cases <- out$prediction
predicted_r <- out$Rt
predicted_deaths <- out$E_deaths

#### Saving each plot individually ####

pdf(file = paste0('results/', indRegion, '_', region,'_cases_', ifrSrc, '.pdf'), width = 8.2, height = 3)
plotInfections(land = region, predicted_cases = predicted_cases, cases = cases, N = N, dates = dates)
dev.off()
pdf(file = paste0('results/', indRegion, '_', region,'_darkfigs_', ifrSrc, '.pdf'), width = 8.2, height = 3)
plotDarkfigures(land = land, predicted_cases = predicted_cases, cases = cases, N = N, dates = dates)
dev.off()
pdf(file = paste0('results/', indRegion, '_', region,'_reproduction_', ifrSrc, '.pdf'), width = 8.2, height = 3)
plotReproduction(predicted_r = predicted_r, N = N, dates = dates)
dev.off()
pdf(file = paste0('results/', indRegion, '_', region,'_deaths_', ifrSrc, '.pdf'), width = 8.2, height = 3)
plotDeaths(predicted_deaths = predicted_deaths, deaths = deaths, dates = dates)
dev.off()

#### Saving all plots in a panel ####

pdf(file = paste0('results/', indRegion, '_', region,'_full_', ifrSrc, '.pdf'), width = 8.2, height = 12)
par(mfrow = c(4,1))
plotDeaths(predicted_deaths = predicted_deaths, deaths = deaths, dates = dates)
plotInfections(land = land, predicted_cases = predicted_cases, cases = cases, N = N, dates = dates)
plotDarkfigures(land = land, predicted_cases = predicted_cases, cases = cases, N = N, dates = dates)
plotReproduction(predicted_r = predicted_r, N = N, dates = dates)
dev.off()
