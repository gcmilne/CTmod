functions {
  real[] si(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    //parameters to fit
    real lambda0 = theta[1];
    int agrps = x_i[1];  //no. age groups
    int K = x_i[2];  //no. state variables
    real r = x_r[1];
    real da = x_r[2];
    real dydt[(agrps*K)];  //define derivative length
    // real d;
    // real propfert[agrps];
    // real deaths[agrps];
    real Na[agrps];
    // real births_age[agrps];
    // real births;
    
    // define states
    real S[agrps]  = y[1:agrps];
    real I[agrps]  = y[(agrps+1):(2*agrps)];
    real Im[agrps] = y[((2*agrps)+1):(3*agrps)];
    
    //total modelled population size
    for(i in 1:agrps){
      Na[i] = y[i] + y[agrps+i] + y[2*agrps+i];
    }
    
    // //age-specific total deaths
    // for(i in 1:agrps){
    //   deaths[i] = d[i] * Na[i]; 
    // }
    // 
    // //age-specific total births
    // for(i in 1:agrps){
    //   births_age[i] = deaths[i] * propfert[i];
    // }
    // 
    // //total births (all ages)
    // births = sum(births_age);
    
    for(i in 1:agrps){
      if(i==1){
        //S
        dydt[i] =  r*Im[i] - lambda0 * S[i] - da*S[i];
        //I
        dydt[agrps+i] = lambda0*(Na[i]-I[i]) - da*I[i];
        //Im
        dydt[2*agrps+i] = -(lambda0 + r + da) * Im[i];
      } else if(i>1){
        //S
        dydt[i] = da*S[i-1] + r*Im[i] - lambda0 * S[i] - da*S[i];
        //I
        dydt[agrps+i] = da*I[i-1] + lambda0*(Na[i]-I[i]) - da*I[i];
        //Im
        dydt[2*agrps+i] = da*Im[i-1] - (lambda0 + r + da) * Im[i];
      } else{
        //S
        dydt[i] = da*S[i-1] + r*Im[i] - lambda0 * S[i];
        //I
        dydt[agrps+i] = da*I[i-1] + lambda0*(Na[i]-I[i]);
        //Im
        dydt[2*agrps+i] = da*Im[i-1] - (lambda0 + r) * Im[i];
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
  // real d;
  // real propfert[agrps];
  real r;
  real da;
  
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
  real x_r[2] = {r, da};   // da, r, propfert, d...
  int x_i[2] = {agrps, K};   //no. age groups, no. state variables
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
