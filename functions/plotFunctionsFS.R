plotInfections <- function(land, predicted_cases, cases, N, dates){
  Sys.setlocale('LC_TIME', 'English') 
  quant_c <- matrix(rep(0, 5*N), nrow = 5)
  for(i in 1:N){
    quant_c[,i] <- quantile(predicted_cases[,i], probs = c(0.025,0.25,0.5,0.75,0.975))
  }
  dates_beg <- max(c(as.Date('2020-02-15'), dates[1]))
  dates_end <- as.Date('2021-01-31')
  dates_seq <- seq.Date(from = dates_beg, to = dates_end, by = 1)
  beg_ind <- which(dates == dates_beg)
  end_ind <- which(dates == dates_end)
  dates_beg <- as.Date('2020-02-15')
  dates_axis <- seq.Date(from = as.Date('2020-02-01'), to = as.Date('2021-02-01'), by = 'month')
  
  offset <- 10
  labs <- match(as.character(dates[beg_ind:end_ind]),as.character(seq.Date(from = as.Date('2019-12-30'), to = max(dates), by = 1)))
  labs <- c(labs[1]-3, labs[1]-2, labs[1]-1, labs, labs[length(labs)]+1, labs[length(labs)]+2, labs[length(labs)]+3)
  impC <- cases[labs]
  meanC <- rep(0, length(dates_seq))
  
  for(i in (4):(length(labs)-3)){
    meanC[i-3] <- mean(impC[((i-3):(i+3))])
  }
  
  par(mar = c(2,6,0,1) + 0.3)
  plot(c(dates_beg, as.Date('2021-02-01')), c(1,1), ann = FALSE, xaxt = 'n', xlab = 'dates', type = 'n', ylim = c(0,max(max(quant_c[,beg_ind:end_ind]),max(cases))), lwd = 0.3, yaxt = 'none', axes = FALSE)
  axis(2, at = pretty(quant_c[5,beg_ind:end_ind]), labels = sprintf('%d',pretty(quant_c[5,beg_ind:end_ind])), cex.axis = 1, lwd = 0.3, las = 2)
  polygon(c(dates[beg_ind:end_ind], rev(dates[beg_ind:end_ind])), c(quant_c[1,beg_ind:end_ind], rev(quant_c[5,beg_ind:end_ind])), col = 'lightblue', border = NA)
  polygon(c(dates[beg_ind:end_ind], rev(dates[beg_ind:end_ind])), c(quant_c[2,beg_ind:end_ind], rev(quant_c[4,beg_ind:end_ind])), col = '#3399FF', border = NA)
  lines(dates[beg_ind:end_ind],quant_c[3,beg_ind:end_ind], lty = 2)
  lines(dates[beg_ind:end_ind],cases[match(as.character(dates[beg_ind:end_ind]),as.character(seq.Date(from = as.Date('2019-12-30'), to = max(dates), by = 1)))], type = 's')
  lines(dates_seq, meanC, type = 'l', col = 'red', lwd = 2)
  labels <- sub('.','',format(as.POSIXct(dates_axis, origin='1970-01-01'), '%d %b %y'))
  axis(1, dates_axis, format(dates_axis, '%b %d'), cex.axis = 0.9, lwd = 0.3, labels = labels, cex = 0.9)
  legend(as.Date('2020-04-01')-7, max(max(quant_c[,beg_ind:end_ind]),max(cases)), c('Posterior median', '95%-CI', '50%-CI', 'Confirmed Cases', '7 day mean Cases'), bty = 'n',
         col = c('black', 'lightblue', '#3399FF', 'black', 'red'), lty = c(2,1,1,1,1), lwd = c(1,2,2,2,2))
  mtext('Infections', side = 2, line = 3.5, cex = 0.8)
}

plotDarkfigures <- function(land, predicted_cases, cases, N, dates){
  Sys.setlocale('LC_TIME', 'English') 
  quant_c <- matrix(rep(0, 5*N), nrow = 5)
  for(i in 1:N){
    quant_c[,i] <- quantile(predicted_cases[,i], probs = c(0.025,0.25,0.5,0.75,0.975))
  }
  dates_beg <- max(c(as.Date('2020-02-15'), dates[1]))
  dates_end <- as.Date('2021-01-31')
  dates_seq <- seq.Date(from = dates_beg, to = dates_end, by = 1)
  beg_ind <- which(dates == dates_beg)
  end_ind <- which(dates == dates_end)
  dates_beg <- as.Date('2020-02-15')
  dates_axis <- seq.Date(from = as.Date('2020-02-01'), to = as.Date('2021-02-01'), by = 'month')
  offset <- 10
  labs <- match(as.character(dates[beg_ind:end_ind]),as.character(seq.Date(from = as.Date('2019-12-30'), to = max(dates), by = 1)))
  labs <- c(labs[1]-3, labs[1]-2, labs[1]-1, labs, labs[length(labs)]+1, labs[length(labs)]+2, labs[length(labs)]+3)
  impC <- cases[labs]
  meanC <- rep(0, length(dates_seq))
  
  quant_cl <- quant_c[,(beg_ind:end_ind)+10]
  darkFigs <- matrix(rep(0, ncol(quant_cl)*nrow(quant_cl)), ncol = ncol(quant_cl))
  
  for(i in (4):(length(labs)-3)){
    meanC[i-3] <- mean(impC[((i-3):(i+3))])
    for(j in 1:nrow(darkFigs)){
      darkFigs[j,i-3] <- quant_cl[j,i-3]/meanC[i-3]
    }
  }
  
  tests <- read.csv('data/rki_tests.CSV', sep = ';', stringsAsFactors = FALSE)
  colnames(tests) <- c('kw', 'sum')
  tests$kw <- seq.Date(from = as.Date('2020-03-02'), by = 'week', length.out = nrow(tests))
  tests <- tests[which(tests$kw <= dates_end),]
  
  par(mar = c(2,6,0,1) + 0.3)
  plot(c(dates_beg, as.Date('2021-02-01')), c(1,1), type = 'n', ylim = c(0,10), ann = FALSE, xaxt = 'n', xlab = 'dates', yaxt = 'none', axes = FALSE)
  polygon(c(dates_seq[which(darkFigs[3,]<8)], rev(dates_seq[which(darkFigs[3,]<8)])), c(darkFigs[1,which(darkFigs[3,]<8)], rev(darkFigs[5,which(darkFigs[3,]<8)])), col = 'lightgrey', border = NA)
  polygon(c(dates_seq[which(darkFigs[3,]<8)], rev(dates_seq[which(darkFigs[3,]<8)])), c(darkFigs[2,which(darkFigs[3,]<8)], rev(darkFigs[4,which(darkFigs[3,]<8)])), col = 'darkgrey', border = NA)
  lines(dates_seq[which(darkFigs[3,]<8)], darkFigs[3,which(darkFigs[3,]<8)], type = 'l', col = 'black')
  lines(dates[c(beg_ind,end_ind)], c(1,1) , lty = 2, col = 'black')
  labels <- sub('.','',format(as.POSIXct(dates_axis, origin='1970-01-01'), '%d %b %y'))
  axis(1, dates_axis, format(dates_axis, '%b %d'), cex.axis = 0.9, lwd = 0.3, labels = labels, cex = 0.9)
  axis(2, pretty(0:10), cex.axis = 1, lwd = 0.3, las = 2)
  legend(as.Date('2020-04-01')-7, 10, c('Posterior median', '95%-CI', '50%-CI'), bty = 'n', pch = c(NA, NA, NA), col = c('black', 'lightgrey', 'darkgrey'), lty = c(1,1,1), lwd = c(2,2,2))
  mtext('Dark Figures', side = 2, line = 3.5, cex = 0.8)
}

plotReproduction <- function(predicted_r, N, dates){
  Sys.setlocale('LC_TIME', 'English') 
  
  quant_r <- matrix(rep(0, 5*N), nrow = 5)
  for(i in 1:N){
    quant_r[,i] <- quantile(predicted_r[,i], probs = c(0.025,0.25,0.5,0.75,0.975))
  }
  
  dates_beg <- max(c(as.Date('2020-02-15'), dates[1]))
  dates_end <- as.Date('2021-01-31')
  dates_seq <- seq.Date(from = dates_beg, to = dates_end, by = 1)
  beg_ind <- which(dates == dates_beg)
  end_ind <- which(dates == dates_end)
  dates_beg <- as.Date('2020-02-15')
  dates_axis <- seq.Date(from = as.Date('2020-02-01'), to = as.Date('2021-02-01'), by = 'month')
  
  par(mar = c(2,6,0,1) + 0.3)
  plot(c(dates_beg, as.Date('2021-02-01')), c(1,1), col = '#FFFFFF', ann = FALSE, xaxt = 'n', xlab = 'dates', type = 's', ylim = c(0,ceiling(min(max(c(max(quant_r),3.5)), 4))+0.2), lwd = 0.3, yaxt = 'none', axes = FALSE)
  labels <- sub('.','',format(as.POSIXct(dates_axis, origin='1970-01-01'), '%d %b %y'))
  axis(1, dates_axis, format(dates_axis, '%b %d'), cex.axis = 0.9, lwd = 0.3, labels = labels, cex = 0.9)
  axis(2, seq(0, ceiling(max(c(max(quant_r),3.5))), 0.5), cex.axis = 1, lwd = 0.3, las = 2)
  polygon(c(dates[beg_ind:end_ind], rev(dates[beg_ind:end_ind])), c(quant_r[1,beg_ind:end_ind], rev(quant_r[5,beg_ind:end_ind])), col = '#C9FFC9', border = NA)
  polygon(c(dates[beg_ind:end_ind], rev(dates[beg_ind:end_ind])), c(quant_r[2,beg_ind:end_ind], rev(quant_r[4,beg_ind:end_ind])), col = '#00FF33', border = NA)
  lines(dates[beg_ind:end_ind],quant_r[3,beg_ind:end_ind], lty = 2)
  lines(dates[beg_ind:end_ind], rep(1, length(dates[beg_ind:end_ind])), col = 'gray', lty = 2)
  legend(as.Date('2020-04-01')-7, ceiling(min(max(c(max(quant_r),3.5)), 4))+0.2, c('Posterior median', '95%-CI', '50%-CI'), 
         bty = 'n', col = c('black', '#C9FFC9', '#00FF33'), lty = c(2,1,1), lwd = c(1,2,2))
  mtext('Reproduction Number', side = 2, line = 3.5, cex = 0.8)
}

plotDeaths <- function(predicted_deaths, deaths, dates){
  Sys.setlocale('LC_TIME', 'English') 
  dates_beg <- max(c(as.Date('2020-02-15'), dates[1]))
  dates_end <- as.Date('2021-01-31')
  dates_seq <- seq.Date(from = dates_beg, to = dates_end, by = 1)
  beg_ind <- which(dates == dates_beg)
  end_ind <- which(dates == dates_end)
  dates_beg <- as.Date('2020-02-15')
  dates_axis <- seq.Date(from = as.Date('2020-02-01'), to = as.Date('2021-02-01'), by = 'month')
  
  quant_d <- matrix(rep(0, 5*N), nrow = 5)
  for(i in 1:N){
    quant_d[,i] <- quantile(predicted_deaths[,i], probs = c(0.025,0.25,0.5,0.75,0.975))
  }
  
  axis_labels <- quant_d[5,beg_ind:end_ind]
  if(max(deaths)>max(axis_labels)) axis_labels <- deaths
  
  par(mar = c(2,6,0,1) + 0.3)
  plot(c(dates_beg, as.Date('2021-02-01')), c(1,1), ann = FALSE, xaxt = 'n', xlab = 'dates', type = 'n', ylim = c(0,max(max(quant_d[,beg_ind:end_ind]),max(deaths))), lwd = 0.3, yaxt = 'none', axes = FALSE)
  polygon(c(dates[beg_ind:end_ind], rev(dates[beg_ind:end_ind])), c(quant_d[1,beg_ind:end_ind], rev(quant_d[5,beg_ind:end_ind])), col = '#FF9966', border = NA)
  polygon(c(dates[beg_ind:end_ind], rev(dates[beg_ind:end_ind])), c(quant_d[2,beg_ind:end_ind], rev(quant_d[4,beg_ind:end_ind])), col = '#FF6633', border = NA)
  lines(dates[beg_ind:end_ind],quant_d[3,beg_ind:end_ind], lty = 2)
  
  lines(dates[beg_ind:end_ind],deaths[match(as.character(dates[beg_ind:end_ind]),as.character(seq.Date(from = as.Date('2019-12-30'), to = max(dates), by = 1)))], type = 's')
  labels <- sub('.','',format(as.POSIXct(dates_axis, origin='1970-01-01'), '%d %b %y'))
  axis(1, dates_axis, format(dates_axis, "%b %d"), cex.axis = 0.9, lwd = 0.3, labels = labels, cex = 0.9)
  axis(2, pretty(axis_labels,7), cex.axis = 1, lwd = 0.3, las = 2)
  legend(as.Date('2020-04-01')-7, max(max(quant_d[,beg_ind:end_ind]),max(deaths)), c('Posterior median', '95%-CI', '50%-CI', 'Reported Deaths'), 
         bty = 'n', col = c('black', '#FF9966', '#FF6633', 'black'), lty = c(2,1,1,1), lwd = c(1,2,2,2))
  mtext('Deaths', side = 2, line = 3.5, cex = 0.8)
}