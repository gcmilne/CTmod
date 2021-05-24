## Define model
age_si = function(time, y, pars) {
  
  # Back-transform parameters
  lambda0  <- exp(pars$log.lambda0)
  lambda1  <- exp(pars$log.lambda1)
  gradient <- exp(pars$log.gradient)
  shape    <- exp(pars$log.shape)
  tdecline <- exp(pars$log.tdecline)
  
  ## set up state variables from input
  S <- y[1:pars$agrps]
  # Infected - either born with congenital disease (seroconversion during pregnancy) or via FoI
  I <- y[(pars$agrps+1):(2*pars$agrps)]
  # Carrier of maternal antibodies at birth b/c of seroconversion in pregnancy but not congenitally diseased
  Im <- y[(2*pars$agrps+1):(3*pars$agrps)]
  
  # set derivatives
  dS <- dI <- dIm <- vector("numeric",length=pars$agrps)
  
  ########################
  ## force of infection ##
  ########################
  # ### set up linear force of infection decrease
  t <- c((pars$burnin - round(tdecline))-1, (pars$burnin+pars$tdiff))     #time period to interpolate
  shape_diff <- c(1, shape)                                               #shape at t[1] and t[2]
  t_shape <- spline(t, shape_diff, xout=seq(t[1], t[2], by=1))$y          #interpolate
  lin_shape <- c( rep(1, length(seq(min(time), (t[1]-1)))), t_shape)
  # plot(time, lin_shape)

  ### set up stepwise foi decrease
  step_shape <- c(rep(1, t[1]), rep(shape, length(t[1]:t[2])-1))

  # create foi array based on whether foi model is constant or age-varying
  if(pars$constant==1){
    foi <- vector("numeric", length=length(time))          #age-constant foi
    
  } else if (pars$constant==0){
    foi <- rep(list(vector("numeric", length=pars$agrps)), length(time)) #age-varying foi
  }

  ## if constant age-foi...
  if(pars$constant==1){
    
    if(pars$stepwise==0){
      foi <- lin_shape * lambda0
      
    }else if (pars$stepwise==1){
      foi <- step_shape * lambda0
    }
    
    ## if age-varying foi...
  } else if (pars$constant==0){
    
    if (pars$stepwise==0){       #linear foi decrease
      
      for(j in 1:length(time)){
        foi[[j]] <- (lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age))) * lin_shape[j]
      }
      
    } else if(pars$stepwise==1){ #stepwise foi decrease
      
      for(j in 1:length(time)){
        foi[[j]] <- (lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age)))*step_shape[j]
      }
      
    }
  }
  
  # calculate change in seroprev and no. seroconversions in pregnancy
  dprev     <- vector("numeric", length=pars$agrps)
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
  
  # modelled population size
  Na <- S+I+Im
  
  ## total deaths
  deaths <- sum((pars$d)*Na)

  ## births distributed among age groups according to fertility
  births_age <-  deaths*pars$propfert
  births <- sum(births_age)
  
  ## move age back 3, 6 or 9 mo to calculate conception distribution for 3 trimesters
  # e.g. c3 = conceived ~9 months ago (more accurately, 7.5 months ago)
  for(i in 1:(pars$agrps-3)){
    c1[i] <- births_age[i+1] 
    c2[i] <- births_age[i+2] 
    c3[i] <- births_age[i+3] 
  }
  
  # prevalence by age (before update)
  pI = (I+Im)/Na
  
  # Adjusting seroprevalence according to equations from Diggle (2011)
  obs_pI <- pI * (pars$se + pars$sp-1) + (1-pars$sp)

  # Calculating seroconversions in pregnancy and cases of congenital disease
  for(i in 1:(pars$agrps-3)){
    if(i==1){
      dprev[i] <- 0
      seroconv1[i] <- 0
      seroconv2[i] <- 0
      seroconv3[i] <- 0
      ct1[i] <- 0
      ct2[i] <- 0
      ct3[i] <- 0
      matAb1[i] <- 0
      matAb2[i] <- 0
      matAb3[i] <- 0
      
    } else {
      dprev[i] <- pI[i]-pI[i-1]                     # change in prevalence (must be positive)
      seroconv1[i] <- dprev[i]*c1[i]                # pregnant women seroconverting in trimester 1
      seroconv2[i] <- dprev[i]*c2[i]                # pregnant women seroconverting in trimester 2
      seroconv3[i] <- dprev[i]*c3[i]                # pregnant women seroconverting in trimester 3
      ct1[i+3] <- seroconv1[i]*pars$mctr[1]         # likelihood of transmission trimester 1
      ct2[i+2] <- seroconv2[i]*pars$mctr[2]         # likelihood of transmission trimester 2
      ct3[i+1] <- seroconv3[i]*pars$mctr[3]         # likelihood of transmission trimester 3
      matAb1[i+3] <- seroconv1[i]*(1-pars$mctr[1])  # maternal Ab trimester 1
      matAb2[i+2] <- seroconv2[i]*(1-pars$mctr[2])  # maternal Ab trimester 2
      matAb3[i+1] <- seroconv3[i]*(1-pars$mctr[3])  # maternal Ab trimester 3
    }
  }
  
  # total number of antibody positive and congenitally diseased births
  matAbt <- sum(matAb1) + sum(matAb2) + sum(matAb3)
  ctt <- sum(ct1) + sum(ct2) + sum(ct3)

  ## select correct model based on assumption of age-invariant foi or age-varying foi
  if(pars$constant==1){
    
    for(j in 1:length(time)){
      for (i in 1:pars$agrps) {
        if (i==1) {
          # susceptible - born seronegative or having lost maternal antibodies, no previous exposure
          dS[i]  <-  (births - matAbt - ctt) + pars$r*Im[i] - foi[j]*S[i] - pars$d[i]*S[i] - pars$da*S[i] 
          # infected - either congenitally or by foi
          dI[i]  <- (ctt + foi[j]*(Na[i]-I[i]) - pars$d[i]*I[i] - pars$da*I[i])
          # maternal antibody positive 
          dIm[i] <- matAbt - (foi[j]+ pars$r+ pars$d[i] + pars$da)*Im[i]
          
        } else if (i<pars$agrps) {
          dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[j]*S[i] - pars$d[i]*S[i] - pars$da*S[i] 
          dI[i]  <- pars$da*I[i-1]  + foi[j]*(Na[i] - I[i]) - pars$d[i] * I[i] - pars$da*I[i]
          dIm[i] <- pars$da*Im[i-1] - (foi[j] + pars$r + pars$d[i] + pars$da) * Im[i] 
          
        } else {
          dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[j]*S[i] - (pars$d[i])*S[i]
          dI[i]  <- pars$da*I[i-1]  + foi[j]*(Na[i]-I[i]) - (pars$d[i])*I[i] 
          dIm[i] <- pars$da*Im[i-1] - (foi[j] + pars$r + pars$d[i])*Im[i] 
          
        }
      }
    }
    
  } else if (pars$constant==0){
    
    for(j in 1:length(time)){
      for (i in 1:pars$agrps) {
        if (i==1) {
          # susceptible - born seronegative or having lost maternal antibodies, no previous exposure
          dS[i]  <-  (births - matAbt - ctt) + pars$r*Im[i] - foi[[j]][i]*S[i] - pars$d[i]*S[i] - pars$da*S[i] 
          # infected - either congenitally or by foi
          dI[i]  <- (ctt + foi[[j]][i]*(Na[i]-I[i]) - pars$d[i]*I[i] - pars$da*I[i])
          # maternal antibody positive 
          dIm[i] <- matAbt - (foi[[j]][i]+ pars$r+ pars$d[i] + pars$da)*Im[i]
          
        } else if (i<pars$agrps) {
          dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[[j]][i]*S[i] - pars$d[i]*S[i] - pars$da*S[i] 
          dI[i]  <- pars$da*I[i-1]  + foi[[j]][i]*(Na[i] - I[i]) - pars$d[i] * I[i] - pars$da*I[i]
          dIm[i] <- pars$da*Im[i-1] - (foi[[j]][i] + pars$r + pars$d[i] + pars$da) * Im[i] 
          
        } else {
          dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[[j]][i]*S[i] - (pars$d[i])*S[i]
          dI[i]  <- pars$da*I[i-1]  + foi[[j]][i]*(Na[i]-I[i]) - (pars$d[i])*I[i] 
          dIm[i] <- pars$da*Im[i-1] - (foi[[j]][i] + pars$r + pars$d[i])*Im[i] 
          
        }
      }
    }
    
  }
  
  ## total susceptible, infected etc
  It <- sum(I)
  Imt <- sum(Im)
  St <- sum(S)
  
  ## recalculate prevalence by age
  pI = (I+Im)/Na
  
  ## recalculate observed prevalence
  obs_pI <- pI * (pars$se + pars$sp-1) + (1-pars$sp)
  
  ## check total pop size
  Nt <- St+It+Imt
  
  ## proportion infected
  pIt <- (It+Imt)/Nt
  
  return(list(y=c(dS, dI, dIm), pI=pI, obs_pI=obs_pI, dprev=dprev, 
              seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
              matAb1=matAb1, matAb2=matAb2, matAb2=matAb2, 
              ct1=ct1, ct2=ct2, ct3=ct3, Na=Na, 
              Nt=Nt, St=St, It=It, Imt=Imt, pIt=pIt, matAbt=matAbt, ctt=ctt))
}
