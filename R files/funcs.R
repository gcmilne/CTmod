###############
## Functions ##
###############

## Function to get age profiles at given time (for 3-month age windows)
if (pars$grps_per_year == 4) {
  
  getit <- function(time) {
    row <- which(abs(sol[,"time"]-time)==min(abs(sol[,"time"]-time)))
    df <- sol[row,-1]
    S         <- df[1:pars$agrps]  
    I         <- df[(pars$agrps+1):(2*pars$agrps)]  
    Im        <- df[(2*pars$agrps+1):(3*pars$agrps)]
    pI        <- df[(3*pars$agrps+1):(4*pars$agrps)]
    dprev     <- df[(4*pars$agrps+1):(5*pars$agrps)]
    seroconv1 <- df[(5*pars$agrps+1):(6*pars$agrps)] 
    seroconv2 <- df[(6*pars$agrps+1):(7*pars$agrps)] 
    seroconv3 <- df[(7*pars$agrps+1):(8*pars$agrps)] 
    matAb1    <- df[(8*pars$agrps+1):(9*pars$agrps)]
    matAb2    <- df[(9*pars$agrps+1):(10*pars$agrps)]
    matAb3    <- df[(10*pars$agrps+1):(11*pars$agrps)]
    ct1       <- df[(11*pars$agrps+1):(12*pars$agrps)]
    ct2       <- df[(12*pars$agrps+1):(13*pars$agrps)]
    ct3       <- df[(13*pars$agrps+1):(14*pars$agrps)]
    Na        <- df[(14*pars$agrps+1):(15*pars$agrps)]
    age       <- pars$age  
    out <- data.frame(a=age, I=I,Im=Im, pI=pI, dprev=dprev, 
                      seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
                      matAb1=matAb1, matAb2=matAb2, matAb3=matAb3, 
                      ct1=ct1, ct2=ct2, ct3=ct3, Na=Na)
    ## remove last age category for aesthetics
    # out <- out[-nrow(out),]
  }
  
} else if (pars$grps_per_year == 12/9) {
  
  ## Function to get age profiles at given time (for 9-month age windows)
  getit <- function(time) {
    row <- which(abs(sol[,"time"]-time)==min(abs(sol[,"time"]-time)))
    df  <- sol[row,-1]
    S        <- df[1:pars$agrps]  
    I        <- df[(pars$agrps+1):(2*pars$agrps)]  
    Im       <- df[(2*pars$agrps+1):(3*pars$agrps)]
    pI       <- df[(3*pars$agrps+1):(4*pars$agrps)]
    dprev    <- df[(4*pars$agrps+1):(5*pars$agrps)]
    seroconv <- df[(5*pars$agrps+1):(6*pars$agrps)] 
    matAb    <- df[(6*pars$agrps+1):(7*pars$agrps)]
    ct       <- df[(7*pars$agrps+1):(8*pars$agrps)]
    Na       <- df[(8*pars$agrps+1):(9*pars$agrps)]
    age      <- pars$age  
    out <- data.frame(a=age, I=I,Im=Im, pI=pI, dprev=dprev, 
                      seroconv=seroconv, matAb=matAb, ct=ct, Na=Na)
    ## remove last age category for aesthetics
    # out <- out[-nrow(out),]
  }
  
}


## Binomial likelihood function ##
loglik <- function (k, n, prev) {
  dbinom(k, n, prev, log = T)
}


## Function to find relevant immunoassay sensitivity and specificity values ## 
# 'fitting_data' seroprevalence data must contain immunoassay method names
# 'assays' immunoassay performance data must contain identical immunoassay method names
find_diagnostic_values <- function(fitting_data, assays) { 
  
  # create df to return
  diagnostic_values <- setNames(data.frame(matrix(nrow=nrow(fitting_data), ncol=2)), c("se", "sp"))
                                
  for(j in 1:nrow(fitting_data)){
    
    ## Select relevant sensitivity & specificity values ##
    if (is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  #if timepoint used 1 immunoassay type
      
      diagnostic_values$se[j] <- assays$se [ which(fitting_data$method[j] == assays$method) ]  #se
      diagnostic_values$sp[j] <- assays$sp [ which(fitting_data$method[j] == assays$method) ]  #sp
      
      
    } else if (!is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  # if timepoint used 2 immunoassays
      
      se_method1 <- assays$se [ which(fitting_data$method [j] == assays$method) ]  #se of first immunoassay
      se_method2 <- assays$se [ which(fitting_data$method2[j] == assays$method) ]  #se of second immunoassay
      
      sp_method1 <- assays$sp [ which(fitting_data$method [j] == assays$method) ]  #sp of first immunoassay
      sp_method2 <- assays$sp [ which(fitting_data$method2[j] == assays$method) ]  #sp of second immunoassay
      
      diagnostic_values$se[j] <- weighted.mean(x = c(se_method1, se_method2), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j])) #weighted mean se
      diagnostic_values$sp[j] <- weighted.mean(x = c(sp_method1, sp_method2), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j])) #weighted mean sp
      
      
    } else if (!is.na(fitting_data$method2[j]) & !is.na(fitting_data$method3[j])) {  # if timepoint used 3 immunoassays
      
      se_method1 <- assays$se [ which(fitting_data$method [j] == assays$method) ]  #se of first immunoassay
      se_method2 <- assays$se [ which(fitting_data$method2[j] == assays$method) ]  #se of second immunoassay
      se_method3 <- assays$se [ which(fitting_data$method3[j] == assays$method) ]  #se of third immunoassay
      
      sp_method1 <- assays$sp [ which(fitting_data$method [j] == assays$method) ]  #sp of first immunoassay
      sp_method2 <- assays$sp [ which(fitting_data$method2[j] == assays$method) ]  #sp of second immunoassay
      sp_method3 <- assays$sp [ which(fitting_data$method3[j] == assays$method) ]  #sp of third immunoassay
      
      diagnostic_values$se[j] <- weighted.mean(x = c(se_method1, se_method2, se_method3), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j], fitting_data$prop_method3[j])) #weighted mean se
      diagnostic_values$sp[j] <- weighted.mean(x = c(sp_method1, sp_method2, sp_method3), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j], fitting_data$prop_method3[j])) #weighted mean sp
      
    }
    
  }
  
  return(diagnostic_values)
}
