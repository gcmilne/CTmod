## Define model
age_si = function(time, y, pars) {
  
  ## set up state variables from input
  # Infected - either born with congenital disease (seroconversion during pregnancy) or via FoI
  S <- y[1:pars$agrps]
  
  # Carrier of maternal antibodies at birth b/c of seroconversion in pregnancy but not congenitally diseased
  I <- y[(pars$agrps+1):(2*pars$agrps)]
  Im <- y[(2*pars$agrps+1):(3*pars$agrps)]
  
  # set derivatives
  dS <- dI <- dIm <- vector("numeric",length=pars$agrps)
  
  ## force of infection
  #burnin <- 250
  #if (time<burnin){
  foi <- (pars$lambda0 + pars$lambda1*(pars$age^2) * (pars$age * exp(-pars$gradient*pars$age)))*pars$shape
  #  } else{
  #  foi <- 0.55*((pars$lambda0 + pars$lambda1*(pars$age^2) * (pars$age * exp(-pars$gradient*pars$age)))*pars$shape)
  #}

  # calculate change in seroprev and no. seroconversions in pregnancy
  dprev     <- vector("numeric", length=pars$agrps)
  seroconv1 <- vector("numeric", length=pars$agrps)
  seroconv2 <- vector("numeric", length=pars$agrps)
  seroconv3 <- vector("numeric", length=pars$agrps)
  matAb1     <- vector("numeric", length=pars$agrps)
  matAb2     <- vector("numeric", length=pars$agrps)
  matAb3     <- vector("numeric", length=pars$agrps)
  c1       <- vector("numeric", length=pars$agrps)
  c2       <- vector("numeric", length=pars$agrps)
  c3       <- vector("numeric", length=pars$agrps)
  ct1       <- vector("numeric", length=pars$agrps)
  ct2       <- vector("numeric", length=pars$agrps)
  ct3       <- vector("numeric", length=pars$agrps)
  
  # modelled population size
  Na <- S+I+Im
  
  ## total deaths
  deaths <- sum(pars$d*Na)
  
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

  for (i in 1:pars$agrps) {
    if (i==1) {
      
      # Susceptible - born seronegative or having lost mat antibodies, no previus exposure
      dS[i]  <-  (births - matAbt - ctt) + pars$r*Im[i] - foi[i]*S[i] - pars$d[i]*S[i] - pars$da*S[i] 
      # Infected - either congenitally or by FoI
      dI[i]  <- (ctt + foi[i]*(Na[i]-I[i]) - pars$d[i]*I[i] - pars$da*I[i])
      # Maternal antibody 
      dIm[i] <- matAbt - (foi[i]+ pars$r+ pars$d[i] + pars$da)*Im[i]
      
    } else if (i<pars$agrps) {
      dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[i]*S[i] - pars$d[i]*S[i] - pars$da*S[i] 
      dI[i]  <- pars$da*I[i-1]  + foi[i]*(Na[i] - I[i]) - pars$d[i] * I[i] - pars$da*I[i]
      dIm[i] <- pars$da*Im[i-1] - (foi[i] + pars$r + pars$d[i] + pars$da) * Im[i] 
      
    } else {
      dS[i]  <- pars$da*S[i-1]  + pars$r*Im[i] - foi[i]*S[i] - (pars$d[i])*S[i]
      dI[i]  <- pars$da*I[i-1]  + foi[i]*(Na[i]-I[i]) - (pars$d[i])*I[i] 
      dIm[i] <- pars$da*Im[i-1] - (foi[i] + pars$r + pars$d[i])*Im[i] 
    }
  }
  
  ## total susceptible, infected etc
  It <- sum(I)
  Imt <- sum(Im)
  St <- sum(S)
  
  ## recalcualte prevalence by age
  pI = (I+Im)/Na
  
  ## check total pop size
  Nt <- St+It+Imt
  
  ## propotion infected
  pIt <- (It+Imt)/Nt
  
  return(list(y=c(dS, dI, dIm), pI=pI, obs_pI=obs_pI, dprev=dprev, 
              seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
              matAb1=matAb1, matAb2=matAb2, matAb2=matAb2, 
              ct1=ct1, ct2=ct2, ct3=ct3, Na=Na, 
              Nt=Nt, St=St, It=It, Imt=Imt, pIt=pIt, matAbt=matAbt, ctt=ctt))
}

