functions {
  real[] si(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    //parameters to fit
    real lambda0 = theta[1];
    int agrps = x_i[1];  //no. age groups
    int K = x_i[2];  //no. state variables
    real dydt[(agrps*K)];  //define derivative length
    
    for(i in 1:agrps){
      //S
      dydt[i] = - lambda0 * y[i];
      //I
      dydt[agrps+i] = lambda0 * y[i];
      //Im
      dydt[2*agrps+i] = lambda0 * y[i];
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
  
  //data to fit
  int<lower=0>cases[agrps];  //no. positive cases
  int<lower=0>n[agrps];  //denominator for each age group
  
  // fixed parameters
  // real r;
  // real mctr[3];
  // real birth_age[agrps];
  
  vector<lower=0>[2] p_lambda0;
  //vector<lower=0>[2] p_seroprev;
  //real p_lambda1;
  //real p_gradient;
  //real p_shape;
  
  //simulation
  real t0;  //starting time
  real ts[200];  //time bins
  int t;  //no. years to run model
  int inference;  //simulate w/o data (inference==0) or with data (inference==1)
  int doprint;
}

transformed data {
  real x_r[0];   // da, r, propfert, d...
  int x_i[2] = {agrps, K};   //no. age groups, no. state variables
  //int x_i[1] = {N};
}

//parameters accepted by the model [that we want values for]
parameters {
  real<lower=0>lambda0;
}

transformed parameters {
  real theta[1];
  // change of format for integrate_ode_rk45
  real init[agrps*K];  // initial values
  real y[t,agrps*K];   // raw ODE outputs
  vector<lower=0, upper=tot_pop>[agrps] comp_S[t];
  vector<lower=0, upper=tot_pop>[agrps] comp_I[t];
  vector<lower=0, upper=tot_pop>[agrps] comp_Im[t];
  vector<lower=0, upper=1>[agrps] comp_pI[t];
  //vector<lower=exp(0)>[agrps] exp_comp_pI[t];  //for likelihood function
  
  // NB by default initial values must be from -2 to 2 on the "unconstrained" scale
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
    ts,      //time bins
    theta,   //parameters
    x_r,     //real data
    x_i,     //integer data
    1.0E-10, //tolerances 
    1.0E-10, //tolerances 
    1.0E3);  //maximum steps
    
    //extract and format ODE results 
    //(1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(i in 1:t) {
      //ok to assume pop strucutre same in S vs. I & Im?
      comp_S[i] = (to_vector(y[i,1:agrps]) + 1.0E-9) * tot_pop;  //total no. (age-stratified) in S
      comp_I[i] = (to_vector(y[i,(agrps+1):(2*agrps)]) + 1.0E-9) * tot_pop;  //total no. (age-stratified) in I
      comp_Im[i] = (to_vector(y[i,(2*agrps+1):(3*agrps)]) + 1.0E-9) * tot_pop;   //total no. (age-stratified) in Im
    }
    
    //compute seroprevalence
    for(i in 1:t) {
      for(j in 1:agrps){
        comp_pI[i,j] = (comp_I[i,j] + comp_Im[i,j]) / comp_S[i,j];
        //exp_comp_pI[i,j] = exp(comp_pI[i,j]);
      }
    }
}

model {
  //priors
  lambda0 ~ lognormal(0, 2); //does log transformation in prior calculation
  
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
    for(i in 1:agrps) {  //for each data point (no. age groups?)
    target += binomial_lpmf(cases[i] | n[i], comp_pI[i]);
    }
  }
}
