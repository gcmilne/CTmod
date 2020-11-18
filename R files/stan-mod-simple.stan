functions{
  vector mod_si(
    real t, vector y,  //y is a vector in new ODE interface 
    real lambda0,  //pass parameters directly into model
    // real lambda1,
    //real gradient, real shape, 
    int agrps, int K, //pass int data directly into model
    // vector age, 
    real r, real da, real[] mctr, real[] d, real[] propfert //pass real data directly into model
    ) {
      
      //define variables to calculate within model 
      real deaths[agrps];
      vector[agrps] Na;
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
      vector[agrps] pI;
      real matAbt;
      real ctt;
      
      //define derivative length
      vector[(agrps*K)] dydt;
      
      // define states
      vector[agrps] S = y[1:agrps];
      vector[agrps] I = y[(agrps+1):(2*agrps)];
      vector[agrps] Im = y[((2*agrps)+1):(3*agrps)];
      
      //define foi
      real foi[agrps];
      foi = rep_array(lambda0, agrps);
      // for(i in 1:agrps){
        //foi[i] = (lambda0 + lambda1*(age[i]^2) * (age[i] * exp(-gradient*age[i])))*shape;
      //   foi[i] = lambda0 + lambda1*(pow(age[i], 2));
      // }
      
      //total modelled population size
      for(i in 1:agrps){
        Na[i] = S[i] + I[i] + Im[i];
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
        pI[i] = (I[i] + Im[i])/Na[i];
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
          dydt[i] = (births - matAbt - ctt) + r*Im[i] - foi[i]*S[i] - d[i]*S[i] - da*S[i];
          //I
          dydt[agrps+i] = ctt + foi[i]*(Na[i]-I[i]) - d[i]*I[i] - da*I[i];
          //Im
          dydt[2*agrps+i] = matAbt - (foi[i] + r + d[i] + da) * Im[i];
          
        } else if(i>1){
          //S
          dydt[i] = da*S[i-1] + r*Im[i] - foi[i]*S[i] - d[i]*S[i] - da*S[i];
          //I
          dydt[agrps+i] = da*I[i-1] + foi[i]*(Na[i]-I[i]) - d[i]*I[i] - da*I[i];
          //Im
          dydt[2*agrps+i] = da*Im[i-1] - (foi[i] + r + d[i] + da) * Im[i];
        } else{
          //S
          dydt[i] = da*S[i-1] + r*Im[i] - foi[i] * S[i] - d[i]*S[i];
          //I
          dydt[agrps+i] = da*I[i-1] + foi[i]*(Na[i]-I[i]) - d[i]*I[i];
          //Im
          dydt[2*agrps+i] = da*Im[i-1] - (foi[i] + r + d[i]) * Im[i];
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
  
  vector[agrps] age;
  real r;
  real da;
  real mctr[3];
  real propfert[agrps];
  real d[agrps];
  
  //data to fit
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
  
  //formatting ode results
  int data_agrps;
  int data_rows[data_agrps*K];
}

// transformed data {
//   real x_r[2+agrps*2] = {r_array};   //r, da, d, propfert
//   int x_i[2] = {agrps, K};   // agrps, K (no. state variables), S_r (no. real arrays), starts_r
// }

//parameters accepted by the model [that we want values for]
parameters {
  real<lower=0, upper=0.2>lambda0;
  //real<lower=0, upper=0.2>lambda1;
  // real<lower=0, upper=0.2>gradient;
  // real<lower=0, upper=0.01>shape;
}

transformed parameters {
  // change of format for integrate_ode_rk45
  vector<lower=0, upper=1>[agrps*K] init;  //initial values
  
  vector[agrps*K] y[t];   //raw ODE outputs
  vector[agrps] comp_S[t];
  vector[agrps] comp_I[t];
  vector[agrps] comp_Im[t];
  vector[agrps] comp_pI[t];
  
  for(i in 1:agrps){
    init[i] = age_prop[i];  //proportion in S0
    init[agrps+i] = 0;      //proportion in I0
    init[2*agrps+i] = 0;    //proportion in Im0
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
    // lambda1, gradient, shape, 
    // age, 
    agrps, K, //pass int data directly into model
    r, da, mctr, d, propfert //pass real data directly into model
    );
    
    //reject simulation if any element of y is NaN
    real na_y[t, K*agrps];
    for(i in 1:t){
      for(j in 1:(3*agrps)){
        na_y[i,j] = is_nan(y[i,j]);
        // reject("this is wrong: y = ", y[i,j]);
      }
    }
    
    //extract and format ODE results
    for(i in 1:t){
      comp_S[i,]  = (to_vector(y[i,1:agrps])) * tot_pop;  //total no. (age-stratified) in S
      comp_I[i,]  = (to_vector(y[i,(agrps+1):(2*agrps)])) * tot_pop;  //total no. (age-stratified) in I
      comp_Im[i,] = (to_vector(y[i,(2*agrps+1):(3*agrps)])) * tot_pop;   //total no. (age-stratified) in Im
    }

    //extract particular age groups specified in the data
    // comp_S  = (to_vector(y[t,data_rows[1:16]])) * tot_pop;  //total no. (age-stratified) in S
    // comp_I  = (to_vector(y[t,data_rows[17:32]])) * tot_pop;  //total no. (age-stratified) in I
    // comp_Im = (to_vector(y[t,data_rows[33:48]])) * tot_pop;   //total no. (age-stratified) in Im
    
    //compute seroprevalence
    for(i in 1:t){
      for(j in 1:agrps){
        comp_pI[i,j] = (comp_I[i,j] + comp_Im[i,j]) / (comp_S[i,j] + comp_I[i,j] + comp_Im[i,j]);
      }
    }
    
}

model {
  lambda0 ~ lognormal(log(.1), .01);
  // lambda1 ~ lognormal(log(.1), .01);
  // gradient ~ lognormal(log(.1), .01);
  // shape ~ lognormal(log(.1), .01);
  
  //debug
  if(doprint==1) {
    // print("na_y: ", na_y[1,]);
    print("lambda0: ", lambda0);
    print("comp_S: ", comp_S[1,]);
    print("comp_I: ", comp_I[1,]);
    print("comp_Im: ", comp_Im[1,]);
    print("comp_pI: ", comp_pI[1,]);
  }
  
  // likelihood 
  if (inference==1) {
    for(i in 1:agrps) { 
    target += binomial_lpmf(cases[i] | n[i], comp_pI[i]);
    }
  }
}

generated quantities {
  vector[agrps] comp_Na[t];
  
  //modelled population size over age and time
  for(i in 1:t){
    for(j in 1:agrps){
      comp_Na[i,j] = comp_S[i,j] + comp_I[i,j] + comp_Im[i,j];
    }
  }
  
}