## Define model
age_si = function(time, y, pars) {
  
  # Back-transform parameters
  lambda0  <- exp(pars$log.lambda0)
  shape    <- exp(pars$log.shape)
  tdecline <- round(exp(pars$log.tdecline), 0)
  
  ## set up state variables from input
  S <- y[ 1:pars$agrps ]
  # Infected - either born with congenital disease (seroconversion during pregnancy) or via FoI
  I <- y[ (pars$agrps + 1) : (2 * pars$agrps) ]
  # Carrier of maternal antibodies at birth b/c of seroconversion in pregnancy but not congenitally diseased
  Im <- y[ (2 * pars$agrps + 1) : (3 * pars$agrps) ]
  
  # set derivatives
  dS <- dI <- dIm <- vector("numeric", length = pars$agrps)
  
  ########################
  ## force of infection ##
  ########################
  
  threshold <- pars$burnin - tdecline
  
  ### Age-constant foi models ###
  if(pars$age_foi == "constant") { 
    
    
    #(1) no temporal decline
    if(pars$temporal_foi == "none") {  
      foi <- lambda0
      
      
      #(2) linear temporal decline
    } else if (pars$temporal_foi == "linear") {  
      
      if (time < threshold) {
        foi <- lambda0 
        
      } else if (time >= threshold) {
        
        yr_rate <- (1 - shape) / ((pars$burnin + pars$tdiff) - threshold)  #define yearly rate of decline
        t_current <- time - threshold  #time difference
        foi <- lambda0 * (1-(yr_rate * t_current)) 
        
        if(foi <= 0){ # to avoid foi < 0 when forecasting
          foi <- 0
        }
      }
      # } else if (time > pars$burnin + pars$tdiff) {  # when forecasting, make foi asymptotic to avoid foi <=0
      #   
      #   yr_rate <- (1 - shape) / ((pars$burnin + pars$tdiff) - threshold)  #define yearly rate of decline
      #   t_current <- time - threshold  #time difference
      #   foi <- (lambda0 * shape)*(exp(-yr_rate*t_current))
      #   
      # }
      
      
      #(3) stepwise temporal decline 
    } else if (pars$temporal_foi == "stepwise") {   
      
      if (time < threshold) {
        foi <- lambda0
        
      } else if (time >= threshold) {
        foi <- lambda0 * shape
      }
      
    }
    
    
    ### Age-doubling foi models ###
  } else if (pars$age_foi == "double") { 
    
    
    #set up foi doubling from ages 20 to 40
    lambda_double <- vector("numeric", length = pars$agrps)  #initiate vector
    lambda_diff <- (2*lambda0 - lambda0)                     #difference in lambda from age group i to i+1
    ages <- (20 * pars$agrps / pars$amax + 1) : (40 * pars$agrps / pars$amax)  #specified age groups
    
    lambda_double <- c(rep(lambda0, (20.5 * pars$agrps / pars$amax) - 1), 
                       seq(lambda0, lambda0 * 2, by = lambda_diff / length(ages)), 
                       rep(lambda0*2, length((40.5 * pars$agrps / pars$amax) : (pars$amax * pars$agrps / pars$amax))))
    
    lambda_double <- lambda_double[-1]
    
    
    #(1) no temporal decline
    if(pars$temporal_foi == "none") {
      foi <- lambda_double
      
      
      #(2) linear temporal decline
    } else if (pars$temporal_foi == "linear") {
      
      if (time < threshold) {
        foi <- lambda_double
        
      } else if (time >= threshold) {
        
        yr_rate <- (1-shape)/((pars$burnin + pars$tdiff) - threshold)  #define yearly rate of decline
        t_current <- time - threshold  #time difference
        foi <- lambda_double * (1-(yr_rate * t_current))
      }
      
      #(3) stepwise temporal decline 
    } else if (pars$temporal_foi == "stepwise") {
      
      if (time < threshold) {
        foi <- lambda_double
        
      } else if (time >= threshold) {
        foi <- lambda_double * shape
      }
      
    }
    
    
    ### Age-halving foi models ###
  } else if (pars$age_foi == "half") {
    
    #set up foi doubling from ages 20 to 40
    lambda_half <- vector("numeric", length=pars$agrps)  #initiate vector
    lambda_diff <- (lambda0 - 0.5 * lambda0)             #difference in lambda from age group i to i+1
    ages <- (20 * pars$agrps/pars$amax + 1) : (40 * pars$agrps / pars$amax)  #specified age groups
    
    lambda_half <- c(rep(lambda0, (20.5 * pars$agrps / pars$amax) - 1), 
                       seq(lambda0, lambda0 * 0.5, by = -lambda_diff / length(ages)), 
                       rep(lambda0 * 0.5, length((40.5 * pars$agrps / pars$amax): (pars$amax * pars$agrps / pars$amax))))
    
    lambda_half <- lambda_half[-1]
    
    
    #(1) no temporal decline
    if(pars$temporal_foi == "none") {
      foi <- lambda_half
      
      
      #(2) linear temporal decline
    } else if (pars$temporal_foi == "linear") {
      
      if (time < threshold) {
        foi <- lambda_half
        
      } else if (time >= threshold) {
        
        yr_rate <- (1 - shape) / ((pars$burnin + pars$tdiff) - threshold)  #define yearly rate of decline
        t_current <- time - threshold  #time difference
        foi <- lambda_half * (1-(yr_rate * t_current))
      }
      
      #(3) stepwise temporal decline 
    } else if (pars$temporal_foi == "stepwise") {
      
      if (time < threshold) {
        foi <- lambda_half
        
      } else if (time >= threshold) {
        foi <- lambda_half * shape
      }
      
    }
    
    
  }
  
  # calculate change in seroprev and no. seroconversions in pregnancy
  dprev     <- vector("numeric", length=pars$agrps)
  
  if(pars$grps_per_year == 4) {
    seroconv1 <- vector("numeric", length=pars$agrps)
    seroconv2 <- vector("numeric", length=pars$agrps)
    seroconv3 <- vector("numeric", length=pars$agrps)
    matAb1    <- vector("numeric", length=pars$agrps)
    matAb2    <- vector("numeric", length=pars$agrps)
    matAb3    <- vector("numeric", length=pars$agrps)
    c1        <- vector("numeric", length=pars$agrps)
    c2        <- vector("numeric", length=pars$agrps)
    c3        <- vector("numeric", length=pars$agrps)
    ct1       <- vector("numeric", length=pars$agrps)
    ct2       <- vector("numeric", length=pars$agrps)
    ct3       <- vector("numeric", length=pars$agrps)
    
  } else if (pars$grps_per_year == 12/9) {
    seroconv <- vector("numeric", length=pars$agrps)
    matAb    <- vector("numeric", length=pars$agrps)
    c_dist   <- vector("numeric", length=pars$agrps)
    ct       <- vector("numeric", length=pars$agrps)
  }
  
  # modelled population size
  Na <- S + I + Im
  
  ## total deaths
  deaths <- sum((pars$d) * Na)
  
  ## total ageing out of max age cat
  byebye <- Na[pars$agrps] * pars$da
  
  ## births distributed among age groups according to fertility (made equal to deaths + ageing beyond max category)
  births_age <- (deaths + byebye) * pars$propfert
  births     <- sum(births_age)
  
  ## Calculate conception distributions 
  if (pars$grps_per_year == 4) { #for 3 trimesters (when age groups = 3 months)
  for (i in 1:(pars$agrps-3)) {
    
    c1[i] <- births_age[i+1] 
    c2[i] <- births_age[i+2] 
    c3[i] <- births_age[i+3]
    
  }
    
  } else if (pars$grps_per_year == 12/9) { #Overall (when age groups = 9 months)
    for (i in 1:(pars$agrps-1)) {
      
      c_dist[i] <- births_age[i+1] 
      
    }
  }
  
  # prevalence by age (before update)
  pI <- (I + Im) / Na
  
  # Adjusting seroprevalence according to equations from Diggle (2011)
  obs_pI <- pI * (pars$se + pars$sp - 1) + (1 - pars$sp)
  
  # Calculating seroconversions in pregnancy and cases of congenital disease
  if (pars$grps_per_year == 4) { # when age groups = 3 months
    
    for (i in 1:(pars$agrps - 3)) {
      
      if (i==1) {
        
        dprev[i]     <- 0
        seroconv1[i] <- 0
        seroconv2[i] <- 0
        seroconv3[i] <- 0
        ct1[i]       <- 0
        ct2[i]       <- 0
        ct3[i]       <- 0
        matAb1[i]    <- 0
        matAb2[i]    <- 0
        matAb3[i]    <- 0
        
      } else {
        
        dprev[i]     <- pI[i] - pI[i-1]                  # change in prevalence (must be positive)
        seroconv1[i] <- dprev[i] * c1[i]                 # pregnant women seroconverting in trimester 1
        seroconv2[i] <- dprev[i] * c2[i]                 # pregnant women seroconverting in trimester 2
        seroconv3[i] <- dprev[i] * c3[i]                 # pregnant women seroconverting in trimester 3
        ct1[i+3]     <- seroconv1[i] * pars$mctr[1]      # likelihood of transmission trimester 1
        ct2[i+2]     <- seroconv2[i] * pars$mctr[2]      # likelihood of transmission trimester 2
        ct3[i+1]     <- seroconv3[i] * pars$mctr[3]      # likelihood of transmission trimester 3
        matAb1[i+3]  <- seroconv1[i] * (1-pars$mctr[1])  # maternal Ab trimester 1
        matAb2[i+2]  <- seroconv2[i] * (1-pars$mctr[2])  # maternal Ab trimester 2
        matAb3[i+1]  <- seroconv3[i] * (1-pars$mctr[3])  # maternal Ab trimester 3
        
      }
      
    }
    
  }  else if (pars$grps_per_year == 12/9) { # when age groups = 9 months
    
    for (i in 1:(pars$agrps - 1)) {
      
      if (i==1) {
        
        dprev[i]     <- 0
        seroconv1[i] <- 0
        c_dist[i]    <- 0
        matAb[i]     <- 0
        
      } else {
        
        dprev[i]     <- pI[i] - pI[i-1]              # change in prevalence (must be positive)
        seroconv1[i] <- dprev[i] * c_dist[i]         # pregnant women seroconverting
        ct[i+1]     <- seroconv1[i] * pars$mctr      # likelihood of transmission
        matAb[i+1]  <- seroconv1[i] * (1-pars$mctr)  # maternal Ab
        
      }
      
    }
    
    
  }

  
  # total number of antibody positive and congenitally diseased births
  if (pars$grps_per_year == 4) {
    matAbt <- sum(matAb1) + sum(matAb2) + sum(matAb3)
    ctt    <- sum(ct1) + sum(ct2) + sum(ct3)
    
  } else if (pars$grps_per_year == 12/9) {
    matAbt <- sum(matAb)
    ctt    <- sum(ct)
  }

  
  
  ## ODEs ##
  
  ## age-constant foi model
  if (pars$age_foi == "constant") { 
    
    for (i in 1:pars$agrps) {
      
      if (i==1) {
        
        # susceptible - born seronegative or having lost maternal antibodies, no previous exposure
        dS[i]  <-  (births - matAbt - ctt) + pars$r*Im[i] - foi*S[i] - pars$d[i]*S[i] - pars$da*S[i]
        
        # infected - either congenitally or by foi
        dI[i]  <- (ctt + foi*(Na[i]-I[i]) - pars$d[i]*I[i] - pars$da*I[i])
        
        # maternal antibody positive 
        dIm[i] <- matAbt - (foi+ pars$r+ pars$d[i] + pars$da)*Im[i]
        

      } else if (i > 1) {

        dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi*S[i] - pars$d[i]*S[i] - pars$da*S[i]
        dI[i]  <- pars$da*I[i-1]  + foi*(Na[i] - I[i]) - pars$d[i] * I[i] - pars$da*I[i]
        dIm[i] <- pars$da*Im[i-1] - (foi + pars$r + pars$d[i] + pars$da) * Im[i]

      }

    }
    
    
    ## age-varying foi models
  } else if (pars$age_foi != "constant") {
    
    for (i in 1:pars$agrps) {
      
      if (i==1) {
        
        # susceptible - born seronegative or having lost maternal antibodies, no previous exposure
        dS[i]  <-  (births - matAbt - ctt) + pars$r*Im[i] - foi[i]*S[i] - pars$d[i]*S[i] - pars$da*S[i]
        
        # infected - either congenitally or by foi
        dI[i]  <- (ctt + foi[i]*(Na[i]-I[i]) - pars$d[i]*I[i] - pars$da*I[i])
        
        # maternal antibody positive
        dIm[i] <- matAbt - (foi[i]+ pars$r+ pars$d[i] + pars$da)*Im[i]
        
      } else if (i <= pars$agrps) {
        
        dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[i]*S[i] - pars$d[i]*S[i] - pars$da*S[i]
        dI[i]  <- pars$da*I[i-1]  + foi[i]*(Na[i] - I[i]) - pars$d[i] * I[i] - pars$da*I[i]
        dIm[i] <- pars$da*Im[i-1] - (foi[i] + pars$r + pars$d[i] + pars$da) * Im[i]
        
      }
    }
    
  }
  
  ## total susceptible, infected etc
  It  <- sum(I)
  Imt <- sum(Im)
  St  <- sum(S)
  
  ## recalculate prevalence by age
  pI  <- (I + Im) / Na
  
  ## recalculate observed prevalence
  obs_pI <- pI * (pars$se + pars$sp - 1) + (1 - pars$sp)
  
  ## check total pop size
  Nt  <- St + It + Imt
  
  ## proportion infected
  pIt <- (It + Imt) / Nt
  
  
  # Sense-check that correct foi model is chosen
  if (pars$troubleshoot == 1) {
    
    if (pars$age_foi != "constant") {
      
      if (time == 1) {
        par(mfrow = c(1,1))
        plot(pars$age, foi, 'l', ylim = c(0, 0.08))
        
      } else if (time > 1 & time < threshold) {
        lines(pars$age, foi)
        
      } else if (time >= threshold) {
        lines(pars$age, foi, lty=2, col = "red")
      }
      
    } else if (pars$age_foi == "constant") {
      
      if (time == 1) {
        par(mfrow = c(1,1))
        plot(pars$age, rep(foi, pars$agrps), 'l', ylim = c(0, 0.08))
        
      } else if (time > 1 & time < threshold) {
        lines(pars$age, rep(foi, pars$agrps))
        
      } else if (time >= threshold) {
        lines(pars$age, rep(foi, pars$agrps), lty=2, col = "red")
      }
      
    }
    
  }
  
  if (pars$grps_per_year == 4) {
    
    return(list(y=c(dS, dI, dIm), pI=pI, obs_pI=obs_pI, dprev=dprev, 
                seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
                matAb1=matAb1, matAb2=matAb2, matAb2=matAb2, 
                ct1=ct1, ct2=ct2, ct3=ct3, Na=Na, 
                Nt=Nt, St=St, It=It, Imt=Imt, pIt=pIt, matAbt=matAbt, ctt=ctt, 
                foi=foi))
    
  } else if (pars$grps_per_year == 12/9) {
    
    return(list(y=c(dS, dI, dIm), pI=pI, obs_pI=obs_pI, dprev=dprev, 
                seroconv=seroconv, matAb=matAb, ct=ct, Na=Na, 
                Nt=Nt, St=St, It=It, Imt=Imt, pIt=pIt, matAbt=matAbt, ctt=ctt, 
                foi=foi))
    
  }

  
}
