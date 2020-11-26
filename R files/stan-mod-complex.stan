functions{
  
  //helper function to find NaNs in a vector
   int count_nans_vec(vector v) {
    int out = 0;
    for(i in 1:cols(v)) {
      if(is_nan(v[i])) {
        out += 1;
      }
    }
    return out;
  }
  
  //function to replace NaNs with 0
  // real zero_nans_real(real r){
  //   int out = 0;
  //   int agrps;
  //   for(i in 1:agrps-3){
  //     if(is_nan(r)==1){
  //       r[i] = 0;
  //     }
  //   }
  //   return out;
  // }
  
  
  //model function
  vector mod_si(
    real t, vector y,  //y is a vector in new ODE interface 
    real lambda0,  //pass parameters directly into model
    // real lambda1,
    // real gradient, real shape,
    int agrps, int K, //pass int data directly into model
    // real[] age, 
    real r, real da, real[] mctr, real[] d, real[] propfert //pass real data directly into model
    ) {
      
      //define variables to calculate within model 
      real deaths[agrps];
      real tot_deaths;
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
      real foi;
      foi = lambda0;
      // for(i in 1:agrps){
      //   foi[i] = (lambda0 + lambda1*(pow(age[i], 2)) * (age[i] * exp(-gradient*age[i])))*shape;
      // }
      
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
      
      // calculating seroconversions in pregnancy and cases of congenital disease
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
          
        } else{
          dprev[i] = pI[i]-pI[i-1];                //change in prevalence (must be positive)
          seroconv1[i] = dprev[i]*c1[i];           //pregnant women seroconverting in trimester 1
          seroconv2[i] = dprev[i]*c2[i];           //pregnant women seroconverting in trimester 2
          seroconv3[i] = dprev[i]*c3[i];           //pregnant women seroconverting in trimester 3
          
          /*   !!  These generate the NaNs    !! */
          //generates 3 NaNs at indices 2:4
          ct1[i+3] = seroconv1[i]*mctr[1];         //likelihood of transmission trimester 1
          //generates 2 NaNs at indices 2:3
          ct2[i+2] = seroconv2[i]*mctr[2];         //likelihood of transmission trimester 2
          //generates 1 NaN at index 2
          ct3[i+1] = seroconv3[i]*mctr[3];         //likelihood of transmission trimester 3
          //generates 3 NaNs at indices 2:4
          matAb1[i+3] = seroconv1[i]*(1-mctr[1]);  //maternal Ab trimester 1
          //generates 2 NaNs at indices 2:3
          matAb2[i+2] = seroconv2[i]*(1-mctr[2]);  //maternal Ab trimester 2
          //generates 1 NaN at index 2
          matAb3[i+1] = seroconv3[i]*(1-mctr[3]);  //maternal Ab trimester 3
        }
      }
      
      //total number of antibody positive and congenitally diseased births
      matAbt = sum(matAb1[5:agrps-3]) + sum(matAb2[4:agrps-3]) + sum(matAb3[3:agrps-3]);
      ctt = sum(ct1[5:agrps-3]) + sum(ct2[4:agrps-3]) + sum(ct3[3:agrps-3]);

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
  
  // real age[agrps];
  real r;
  real da;
  real mctr[3];
  real propfert[agrps];
  real d[agrps];
  
  //formatting ode results
  int data_agrps;
  int data_rows[data_agrps*K];
  
  //data to fit
  int<lower=0>cases[data_agrps];  //no. positive cases
  int<lower=0>n[data_agrps];  //denominator for each age group
  
  //simulation
  real t0;  //starting time
  int t;  //no. years to run model
  real ts[t];  //time bins
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  int inference;  //simulate w/o data (inference==0) or with data (inference==1)
  int doprint;
<<<<<<< HEAD
=======
  
  //formatting ode results
  // int data_agrps;
  // int data_rows[data_agrps*K];
>>>>>>> 8f66a06d9f41818f9e74a208cb90782eb2274b5a
}

parameters {
  real<lower=0, upper=0.10>lambda0;
  // real<lower=0, upper=0.2>lambda1;
  // real<lower=0, upper=0.2>gradient;
  // real<lower=0, upper=0.01>shape;
}

transformed parameters {
  // change of format for integrate_ode_rk45
  vector<lower=0, upper=1>[agrps*K] init;  //initial values
  
<<<<<<< HEAD
  vector[agrps*K] y[t];
  vector[data_agrps] comp_S;
  vector[data_agrps] comp_I;
  vector[data_agrps] comp_Im;
  vector[data_agrps] comp_Na;
  // vector<lower=0, upper=1>[agrps] comp_pI;
  vector[data_agrps] comp_pI;
=======
  vector[agrps*K] y[t];   //raw ODE outputs
  vector[agrps] comp_S[t];
  vector[agrps] comp_I[t];
  vector[agrps] comp_Im[t];
  vector[agrps] comp_Na;
  // vector<lower=0, upper=1>[agrps] comp_pI;
  vector[agrps] comp_pI;
>>>>>>> 8f66a06d9f41818f9e74a208cb90782eb2274b5a
  
  for(i in 1:(agrps)){
    init[i]         = age_prop[i];
    init[agrps+i]   = 0;
    init[2*agrps+i] = 0;
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
    agrps, K, //pass int data directly into model
    // age,
    r, da, mctr, d, propfert //pass real data directly into model
    );
    
    //reject simulation if any element of y is NaN
    // real na_y[t, K*agrps];
    // for(i in 1:t){
    //   for(j in 1:(3*agrps)){
    //     na_y[i,j] = is_nan(y[i,j]);
    //     reject("this is wrong: y = ", y[i,j]);
    //   }
    // }
    
    //extract and format ODE results
<<<<<<< HEAD
    // for(i in 1:t){
    //   for(j in 1:agrps){
    //     comp_S[i,j]  = (to_vector(y[i,j])) * tot_pop;  //total no. (age-stratified) in S
    //     comp_I[i,j]  = (to_vector(y[i,agrps+j])) * tot_pop;  //total no. (age-stratified) in I
    //     comp_Im[i,j] = (to_vector(y[i,2*agrps+j])) * tot_pop;   //total no. (age-stratified) in Im
    //   }
    // }
    // 
=======
    for(i in 1:t){
      comp_S[i,]  = (to_vector(y[i,1:agrps])) * tot_pop;  //total no. (age-stratified) in S
      comp_I[i,]  = (to_vector(y[i,(agrps+1):(2*agrps)])) * tot_pop;  //total no. (age-stratified) in I
      comp_Im[i,] = (to_vector(y[i,(2*agrps+1):(3*agrps)])) * tot_pop;   //total no. (age-stratified) in Im
    }

    //extract particular age groups specified in the data
    comp_S  = (to_vector(y[t,data_rows[1:data_agrps]])) * tot_pop;                    //total no. (age-stratified) in S
    comp_I  = (to_vector(y[t,data_rows[(data_agrps+1):(2*data_agrps)]])) * tot_pop;   //total no. (age-stratified) in I
    comp_Im = (to_vector(y[t,data_rows[(2*data_agrps+1):(3*data_agrps)]])) * tot_pop; //total no. (age-stratified) in Im
    
    //compute seroprevalence
    comp_pI = (comp_I + comp_Im) ./ (comp_S + comp_I + comp_Im);
    
    // //compute Na
    comp_Na = comp_Im + comp_I + comp_S;
    
}

model {
  lambda0 ~ lognormal(log(.045), .3);
  // lambda1 ~ lognormal(log(.1), .01);
  // gradient ~ lognormal(log(.1), .01);
  // shape ~ lognormal(log(.1), .01);
  
  //debug
  if(doprint==1) {
    // print("na_y: ", na_y[1,]);
    print("lambda0: ", lambda0);
    print("comp_S: ", comp_S);          //make sure these look reasonable
    print("comp_I: ", comp_I);
    print("comp_Im: ", comp_Im);
    print("comp_Na: ", comp_Na);
    print("comp_pI: ", comp_pI);
  }
  
  // likelihood 
  if (inference==1) {
<<<<<<< HEAD
    for(j in 1:data_agrps) { 
    target += binomial_lpmf(cases[j] | n[j], comp_pI[j]);
    }
  }
}

// generated quantities {
//   vector[agrps] comp_Na[t];
//   
//   //modelled population size over age and time
//   for(i in 1:t){
//     for(j in 1:agrps){
//       comp_Na[i,j] = comp_S[i,j] + comp_I[i,j] + comp_Im[i,j];
//     }
//   }
//   
// }
