functions{
  
  //helper function to find NaNs
  int count_nans_vec(vector v) {
    int out = 0;
    for(i in 1:cols(v)) {
      if(is_nan(v[i])) {
        out += 1;
      }
    }
    return out;
  }
  
  //model function
  vector mod_si(
    real t, vector y,  //y is a vector in new ODE interface 
    real lambda0,  //pass parameters directly into model
    int agrps, int K, //pass int data directly into model
    real r, real da, real mean_mctr, real[] d, real[] propfert //pass real data directly into model
    ) {
      
      //define variables to calculate within model 
      real deaths[agrps];
      real tot_deaths;
      vector[agrps] Na;
      real births_age[agrps];
      real births;
      real dprev[agrps];
      real seroconv1[agrps];
      real c1[agrps];
      real ct1[agrps];
      real matAb1[agrps];
      vector[agrps] pI;
      real matAbt;
      real ctt;

      //define derivative length
      vector[agrps*K] dydt;
      
      // define states
      vector[agrps] S = y[1:agrps];
      vector[agrps] I = y[(agrps+1):(2*agrps)];
      vector[agrps] Im = y[((2*agrps)+1):(3*agrps)];
      
      //define foi
      real foi = lambda0;
      
      //total modelled population size
      for(i in 1:agrps){
        Na[i] = S[i] + I[i] + Im[i];
      }
      
      //age-specific deaths
      for(i in 1:agrps){
        deaths[i] = d[i] * Na[i];
      }
      
      //total deaths
      tot_deaths = sum(deaths);
      
      //age-specific total births
      for(i in 1:agrps){
        births_age[i] = tot_deaths * propfert[i];
      }
      
      //total births (all ages)
      births = sum(births_age);
      
      // conception distribution (assuming 1 yr pregnancy)
      for(i in 1:(agrps-1)){
        c1[i] = births_age[i+1];
      }

      //seroprevalence
      for(i in 1:agrps){
        pI[i] = (I[i] + Im[i])/Na[i];
      }

      // calculating seroconversions in pregnancy and cases of congenital disease
      for(i in 1:(agrps-1)){
        if(i==1){
          dprev[i]     = 0;
          seroconv1[i] = 0;
          ct1[i]       = 0;
          matAb1[i]    = 0;

        } else{
          dprev[i]     = pI[i]-pI[i-1];            //change in prevalence (must be positive)
          seroconv1[i] = dprev[i]*c1[i];           //pregnant women seroconverting
          ct1[i+1] = seroconv1[i]*mean_mctr;         //likelihood of transmission trimester 1
          matAb1[i+1] = seroconv1[i]*(1-mean_mctr);  //maternal Ab trimester 1
        }
      }

      //total number of antibody positive and congenitally diseased births
      //remove NaNs from summed variables
      matAbt = sum(matAb1[3:agrps-1]);
      ctt    = sum(ct1[3:agrps-1]);

      //model ODEs
      for(i in 1:agrps){
        if(i==1){
          dydt[i] = (births - matAbt - ctt) + r*Im[i] - foi*S[i] - d[i]*S[i] - da*S[i];
          dydt[agrps+i] = ctt + foi*(Na[i]-I[i]) - d[i]*I[i] - da*I[i];
          dydt[2*agrps+i] = matAbt - (foi + r + d[i] + da) * Im[i];
          
        } else if(i>1){
          dydt[i] = da*S[i-1] + r*Im[i] - foi*S[i] - d[i]*S[i] - da*S[i];
          dydt[agrps+i] = da*I[i-1] + foi*(Na[i]-I[i]) - d[i]*I[i] - da*I[i];
          dydt[2*agrps+i] = da*Im[i-1] - (foi + r + d[i] + da) * Im[i];
          
        } else{
          dydt[i] = da*S[i-1] + r*Im[i] - foi * S[i] - d[i]*S[i];
          dydt[agrps+i] = da*I[i-1] + foi*(Na[i]-I[i]) - d[i]*I[i];
          dydt[2*agrps+i] = da*Im[i-1] - (foi + r + d[i]) * Im[i];
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
  
  real r;
  real da;
  real mean_mctr;
  real propfert[agrps];
  real d[agrps];
  
  //formatting ode results
  // int data_agrps;
  // int data_rows[data_agrps*K];
  
  //data to fit
  // int<lower=0>cases[data_agrps];  //no. positive cases
  // int<lower=0>n[data_agrps];  //denominator for each age group
  int<lower=0>cases[agrps];  //no. positive cases
  int<lower=0>n[agrps];  //denominator for each age group

  //simulation
  real t0;  //starting time
  int t;  //no. years to run model
  real ts[t];  //time bins
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  int inference;  //simulate w/o data (inference==0) or with data (inference==1)
  int doprint;
}

parameters {
  real<lower=0, upper=0.10>lambda0;
}

transformed parameters {
  // change of format for integrate_ode_rk45
  vector<lower=0, upper=1>[agrps*K] init;  //initial values
  
  vector[agrps*K] y[t];   //raw ODE outputs
  vector[agrps] comp_S;
  vector[agrps] comp_I;
  vector[agrps] comp_Im;
  vector[agrps] comp_Na;
  vector[agrps] comp_pI;
  // vector[data_agrps] comp_S;
  // vector[data_agrps] comp_I;
  // vector[data_agrps] comp_Im;
  // vector[data_agrps] comp_Na;
  // vector[data_agrps] comp_pI;
  
  for(i in 1:agrps){
    init[i]          = age_prop[i];  //proportion in S0
    init[agrps+i]    = 0;            //proportion in I0
    init[2*agrps+i]  = 0;            //proportion in Im0
  }
  
  //run solver
  y = ode_rk45_tol(
    mod_si,  //model function
    init,    //vector initial values
    t0,      //real initial time
    ts,      //real[] times
    rel_tol,
    abs_tol,
    max_num_steps,
    lambda0,  //pass parameters directly into model
    agrps, K, //pass int data directly into model
    r, da, mean_mctr, d, propfert //pass real data directly into model
    );
    
    //extract and format ODE results
    comp_S  = (to_vector(y[t,1:agrps])) * tot_pop;  
    comp_I  = (to_vector(y[t,(agrps+1):(2*agrps)])) * tot_pop; 
    comp_Im = (to_vector(y[t,(2*agrps+1):(3*agrps)])) * tot_pop; 
    
    //extract and format ODE results given age groups specified in the data
    // comp_S  = (to_vector(y[t,data_rows[1:data_agrps]])) * tot_pop;  
    // comp_I  = (to_vector(y[t,data_rows[(data_agrps+1):(2*data_agrps)]])) * tot_pop; 
    // comp_Im = (to_vector(y[t,data_rows[(2*data_agrps+1):(3*data_agrps)]])) * tot_pop; 
    
    //compute seroprevalence
    comp_pI = (comp_I + comp_Im) ./ (comp_S + comp_I + comp_Im);
    
    // //compute Na
    comp_Na = comp_Im + comp_I + comp_S;
}

model {
  lambda0 ~ lognormal(log(.1), .01);
  
  //debug
  if(doprint==1) {
    print("lambda0: ", lambda0);
    print("comp_S: ", comp_S);
    print("comp_I: ", comp_I);
    print("comp_Im: ", comp_Im);
    print("comp_Na: ", comp_Na);
    print("comp_pI: ", comp_pI);
  }
  
  // likelihood 
  if (inference==1) {
    for(j in 1:agrps) { 
    target += binomial_lpmf(cases[j] | n[j], comp_pI[j]);
    }
  }
}
