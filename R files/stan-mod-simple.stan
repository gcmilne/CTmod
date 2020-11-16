functions {
  real[] si(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    //define parameters
    real lambda0 = theta[1];
    
    //define integer data
    int agrps = x_i[1];  //no. age groups
    int K = x_i[2];  //no. state variables
    int da;
    
    //define real data
    int S_r;  
    int starts_r[S_r + 1]; //start indices
    real r[1] = x_r[starts_r[1]:(starts_r[2] - 1)];              //1st array, r
    // real da[1] = x_r[starts_r[2]:(starts_r[3] - 1)];            //2nd array, da
    real mctr[3] = x_r[starts_r[3]:(starts_r[4] - 1)];          //3rd array, mctr
    real d[agrps] = x_r[starts_r[4]:(starts_r[5] - 1)];         //4th array, d
    real propfert[agrps] = x_r[starts_r[5]:(starts_r[6] - 1)];  //5th array, propfert
    
    // real r_mod[agrps] = rep_array(r, agrps);
    // int da_mod[agrps] = rep_array(da, agrps);
    
    //define variables to calculate within model 
    real deaths[agrps];
    real Na[agrps];
    real births_age[agrps];
    real births;
    real dprev[agrps];
    real seroconv1[agrps];
    real seroconv2[agrps];
    real seroconv3[agrps];
    real c1[agrps];
    real c2[agrps];
    real c3[agrps];
    real ct1[agrps];
    real ct2[agrps];
    real ct3[agrps];
    real matAb1[agrps];
    real matAb2[agrps];
    real matAb3[agrps];
    real pI[agrps];
    real matAbt;
    real ctt;
    
    //define derivative length
    real dydt[(agrps*K)];
        
    // define states
    real S[agrps]  = y[1:agrps];
    real I[agrps]  = y[(agrps+1):(2*agrps)];
    real Im[agrps] = y[((2*agrps)+1):(3*agrps)];
    
    //total modelled population size
    for(i in 1:agrps){
      Na[i] = y[i] + y[agrps+i] + y[2*agrps+i];
    }
    
    //age-specific total deaths
    for(i in 1:agrps){
      deaths[i] = d[i] * Na[i];
    }

    //age-specific total births
    for(i in 1:agrps){
      births_age[i] = deaths[i] * propfert[i];
    }

    //total births (all ages)
    births = sum(births_age);
    
    // conception distribution
    //move age back 3, 6 or 9 mo to calculate conception distribution for 3 trimesters
    //e.g. c3 = conceived ~9 months ago (more accurately, 7.5 months ago)
    for(i in 1:(agrps-3)){
      c1[i] = births_age[i+1];
      c2[i] = births_age[i+2];
      c3[i] = births_age[i+3];
    }

    //seroprevalence
    for(i in 1:agrps){
      pI[i] = (y[agrps+i] + y[2*agrps+i])/Na[i];
    }

    //calculating seroconversions in pregnancy and cases of congenital disease
    for(i in 1:(agrps-3)){
      if(i==1){
        dprev[i] = 0;
        seroconv1[i] = 0;
        seroconv2[i] = 0;
        seroconv3[i] = 0;
        ct1[i] = 0;
        ct2[i] = 0;
        ct3[i] = 0;
        matAb1[i] = 0;
        matAb2[i] = 0;
        matAb3[i] = 0;

      } else {
        dprev[i] = pI[i]-pI[i-1];                //change in prevalence (must be positive)
        seroconv1[i] = dprev[i]*c1[i];           //pregnant women seroconverting in trimester 1
        seroconv2[i] = dprev[i]*c2[i];           //pregnant women seroconverting in trimester 2
        seroconv3[i] = dprev[i]*c3[i];           //pregnant women seroconverting in trimester 3
        ct1[i+3] = seroconv1[i]*mctr[1];         //likelihood of transmission trimester 1
        ct2[i+2] = seroconv2[i]*mctr[2];         //likelihood of transmission trimester 2
        ct3[i+1] = seroconv3[i]*mctr[3];         //likelihood of transmission trimester 3
        matAb1[i+3] = seroconv1[i]*(1-mctr[1]);  //maternal Ab trimester 1
        matAb2[i+2] = seroconv2[i]*(1-mctr[2]);  //maternal Ab trimester 2
        matAb3[i+1] = seroconv3[i]*(1-mctr[3]);  //maternal Ab trimester 3
      }
    }

    //total number of antibody positive and congenitally diseased births
    matAbt = sum(matAb1) + sum(matAb2) + sum(matAb3);
    ctt = sum(ct1) + sum(ct2) + sum(ct3);
    
    //model ODEs
    for(i in 1:agrps){
      if(i==1){
        //S
        dydt[i] = (births - matAbt - ctt) + r*Im[i] - lambda0*S[i] - d[i]*S[i] - da*S[i];
        //I
        dydt[agrps+i] = ctt + lambda0*(Na[i]-I[i]) - d[i]*I[i] - da*I[i];
        //Im
        dydt[2*agrps+i] = matAbt - (lambda0 + r + d[i] + da) * Im[i];
        
      } else if(i>1){
        //S
        dydt[i] = da*S[i-1] + r*Im[i] - lambda0*S[i] - d[i]*S[i] - da*S[i];
        //I
        dydt[agrps+i] = da*I[i-1] + lambda0*(Na[i]-I[i]) - d[i]*I[i] - da*I[i];
        //Im
        dydt[2*agrps+i] = da*Im[i-1] - (lambda0 + r + d[i] + da) * Im[i];
      } else{
        //S
        dydt[i] = da*S[i-1] + r*Im[i] - lambda0 * S[i] - d[i]*S[i];
        //I
        dydt[agrps+i] = da*I[i-1] + lambda0*(Na[i]-I[i]) - d[i]*I[i];
        //Im
        dydt[2*agrps+i] = da*Im[i-1] - (lambda0 + r + d[i]) * Im[i];
      }
    }
    return(dydt);
  }
}

data {
  //for model calculations
  int K; //no. state variables
  int agrps;  //no. age groups
  vector[agrps] age_prop;  //proportion of population in each age group
  int tot_pop;  //total population size
  
  /* building x_r array */
  //no. real arrays: {r, da, mctr, d, propfert} = 5
  real r_array[1+1+3+agrps*2];
  
  
  // /* building x_i array */
  // real i_array[2+(4*2)];  //array to hold all values for x_i
  // 
  // //no. integer arrays: {agrps, K, S_r, starts_r, S_i, starts_i} = 6
  // int S_i;  
  // 
  // int starts_i[S_i + 1]; //start indices of S_i arrays (with the last index being N + 1)
  
  //data to fit
  int<lower=0>cases[agrps];  //no. positive cases
  int<lower=0>n[agrps];  //denominator for each age group

  //simulation
  real t0;  //starting time
  int t;  //no. years to run model
  real ts[t];  //time bins
  int inference;  //simulate w/o data (inference==0) or with data (inference==1)
  int doprint;
}

transformed data {
  real x_r[2+agrps*2] = {r_array};   //r, da, d, propfert
  int x_i[2] = {agrps, K};   // agrps, K (no. state variables), S_r (no. real arrays), starts_r
}

//parameters accepted by the model [that we want values for]
parameters {
  real<lower=0, upper=0.2>lambda0;
}

transformed parameters {
  real theta[1];
  // change of format for integrate_ode_rk45
  real init[agrps*K];  // initial values
  real y[t,agrps*K];   // raw ODE outputs
  vector[agrps] comp_S;
  vector[agrps] comp_I;
  vector[agrps] comp_Im;
  vector[agrps] comp_pI;
  
  theta[1] = lambda0;
  for(i in 1:agrps){
    init[i] = age_prop[i];  //proportion in S0
    init[agrps+i] = 0;      //proportion in I0
    init[2*agrps+i] = 0;    //proportion in Im0
  }
  
  //run solver
  y = integrate_ode_rk45(
    si,      //model fucntion
    init,    //initial values
    t0,      //initial time
    ts,      //evaluation dates
    theta,   //parameters
    x_r,     //real data
    x_i,     //integer data
    1.0E-10, //tolerances 
    1.0E-10, //tolerances 
    1.0E3);  //maximum steps
    
    //extract and format ODE results 
    comp_S = (to_vector(y[t,1:agrps])) * tot_pop;  //total no. (age-stratified) in S
    comp_I = (to_vector(y[t,(agrps+1):(2*agrps)])) * tot_pop;  //total no. (age-stratified) in I
    comp_Im = (to_vector(y[t,(2*agrps+1):(3*agrps)])) * tot_pop;   //total no. (age-stratified) in Im
    
    //compute seroprevalence
    for(j in 1:agrps){
      comp_pI[j] = (comp_I[j] + comp_Im[j]) / (comp_I[j] + comp_Im[j] + comp_S[j]);
    }
}

model {
  lambda0 ~ lognormal(log(.1), .01);
  
  //debug
  if(doprint==1) {
    print("lambda0: ", lambda0);
    print("comp_S: ", comp_S);
    print("comp_I: ", comp_I);
    print("comp_Im: ", comp_Im);
    print("comp_pI: ", comp_pI);
  }
  
  // likelihood 
  if (inference==1) {
    for(i in 1:agrps) { 
    target += binomial_lpmf(cases[i] | n[i], comp_pI[i]);
    }
  }
}
