functions {
  real[] si(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    //parameters to fit
    real lambda0 = theta[1];
    
    int agrps = x_i[1];  //no. age groups
    int K = x_i[2];  //no. state variables
    real r = x_r[1];
    real da = x_r[2];
    real dydt[(agrps*K)];  //define derivative length

    //model running variables
    // real dprev[agrps];
    // real seroconv1[agrps];
    // real seroconv2[agrps];
    // real seroconv3[agrps];
    // real c1[agrps];
    // real c2[agrps];
    // real c3[agrps];
    // real ct1[agrps];
    // real ct2[agrps];
    // real ct3[agrps];
    // real matAb1[agrps];
    // real matAb2[agrps];
    // real matAb3[agrps];
    // real pI[agrps];

    //foi
    // real foi[agrps] = rep_array(lambda0, agrps);  //array of length agrps, with all values set to lambda0
    real foi = lambda0;  //array of length agrps, with all values set to lambda0

    //define states
    // real S[agrps]  = y[1:agrps];
    // real I[agrps]  = y[(agrps+1):(2*agrps)];
    // real Im[agrps] = y[((2*agrps)+1):(3*agrps)];
    
    
    //define births, deaths & congenital infections
    // real mctr[3];
    real d[agrps];
    real propfert[agrps];
    real deaths[agrps];
    // real Na[agrps];
    real births_age[agrps];
    real births;
    real sum_deaths;
    // real matAbt;
    // real ctt;

    //total modelled population size
    // for(i in 1:agrps){
    //   Na[i] = y[i] + y[agrps+i] + y[2*agrps+i];
    // }
    
    //age-specific total deaths
    for(i in 1:agrps){
      // deaths[i] = d[i] * Na[i]; 
      deaths[i] = d[i] * (y[i] + y[agrps+i] + y[2*agrps+i]);
    }
    
    //sum deaths
    sum_deaths = sum(deaths);
    
    //age-specific total births
    for(i in 1:agrps){
      births_age[i] = sum_deaths * propfert[i];
    }
    
    //total births (all ages)
    births = sum(births_age);
    
    //conception distribution
    //move age back 3, 6 or 9 mo to calculate conception distribution for 3 trimesters
    //e.g. c3 = conceived ~9 months ago (more accurately, 7.5 months ago)
    // for(i in 1:(agrps-3)){
    //   c1[i] = births_age[i+1];
    //   c2[i] = births_age[i+2];
    //   c3[i] = births_age[i+3];
    // }
    
    //seroprevalence
    // for(i in 1:agrps){
    //   pI[i] = (y[agrps+i] + y[2*agrps+i])/Na[i];
    // }

    // //calculating seroconversions in pregnancy and cases of congenital disease
    // for(i in 1:(agrps-3)){
    //   if(i==1){
    //     dprev[i] = 0;
    //     seroconv1[i] = 0;
    //     seroconv2[i] = 0;
    //     seroconv3[i] = 0;
    //     ct1[i] = 0;
    //     ct2[i] = 0;
    //     ct3[i] = 0;
    //     matAb1[i] = 0;
    //     matAb2[i] = 0;
    //     matAb3[i] = 0;
    // 
    //   } else {
    //     dprev[i] = pI[i]-pI[i-1];                //change in prevalence (must be positive)
    //     seroconv1[i] = dprev[i]*c1[i];           //pregnant women seroconverting in trimester 1
    //     seroconv2[i] = dprev[i]*c2[i];           //pregnant women seroconverting in trimester 2
    //     seroconv3[i] = dprev[i]*c3[i];           //pregnant women seroconverting in trimester 3
    //     ct1[i+3] = seroconv1[i]*mctr[1];         //likelihood of transmission trimester 1
    //     ct2[i+2] = seroconv2[i]*mctr[2];         //likelihood of transmission trimester 2
    //     ct3[i+1] = seroconv3[i]*mctr[3];         //likelihood of transmission trimester 3
    //     matAb1[i+3] = seroconv1[i]*(1-mctr[1]);  //maternal Ab trimester 1
    //     matAb2[i+2] = seroconv2[i]*(1-mctr[2]);  //maternal Ab trimester 2
    //     matAb3[i+1] = seroconv3[i]*(1-mctr[3]);  //maternal Ab trimester 3
    //   }
    // }

    //total number of antibody positive and congenitally diseased births
    // matAbt = sum(matAb1) + sum(matAb2) + sum(matAb3);
    // ctt = sum(ct1) + sum(ct2) + sum(ct3);
    
    //model
    for(i in 1:agrps){
        //S[i] or y[i]
        //S[i-1] or y[i-1]
        //I[i], or y[agrps+i]
        //I[i-1] or y[(agrps+i-1]
        //Im[i] or y[(2*agrps+i]
      if(i==1){
        //S
        dydt[i] = (births - matAbt - ctt) + r*y[2*agrps+i] - foi[i]*y[i] - d[i]*y[i] - da*y[i];
        // dydt[i] = (births) + r*y[2*agrps+i] - foi*y[i] - d[i]*y[i] - da*y[i];
        //I - infected - either congenitally or by FoI
        // dydt[agrps+i] = ctt + foi[i]*(Na[i]-y[agrps+i]) - d[i]*y[agrps+i] - da*y[agrps+i];
        dydt[agrps+i] = foi*((y[i] + y[agrps+i] + y[2*agrps+i])-y[agrps+i]) - d[i]*y[agrps+i] - da*y[agrps+i];
        //Im - maternal antibody
        // dydt[2*agrps+i] = matAbt - (foi[i]+ r + d[i] + da) * y[2*agrps+i];
        dydt[2*agrps+i] = - (foi+ r + d[i] + da) * y[2*agrps+i];

      }else if(i <agrps){
        //S
        dydt[i] = da*y[i-1]  + r*y[2*agrps+i] - foi*y[i] - d[i]*y[i] - da*y[i];
        //I
        dydt[agrps+i] = da*y[agrps+i-1]  + foi*((y[i] + y[agrps+i] + y[2*agrps+i]) - y[agrps+i]) - d[i]*y[agrps+i] - da*y[agrps+i];
        //Im
        dydt[2*agrps+i] = da*y[2*agrps+i-1] - (foi + r + d[i] + da) * y[2*agrps+i];
        
      } else{
        //S
        dydt[i] =  da*y[i-1]  + r*y[2*agrps+i] - foi*y[i] - d[i]*y[i];
        //I
        dydt[agrps+i] =  da*y[agrps+i-1]  + foi*((y[i] + y[agrps+i] + y[2*agrps+i])-y[agrps+i]) - d[i]*y[agrps+i];
        //Im
        dydt[2*agrps+i] = da*y[2*agrps+i-1] - (foi + r + d[i]) * y[2*agrps+i];
      }
    }
    return dydt;
  }
}

data {
  //for model function
  int K; //no. state variables
  int agrps;  //no. age groups in model
  // real mctr[3];
  real d[agrps];
  // real propfert[agrps];
  real r;
  real da;
  
  //for formatting ODE results
  int data_agrps;  //no. age groups in data
  vector[agrps] age_prop;  //proportion of population in each age group
  int tot_pop;  //total population size
  int data_rows[K*data_agrps];  //rows of df in which data exist
  
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
  real x_r[2] = {r, da};   
  int x_i[2] = {agrps, K};
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
  vector[data_agrps] comp_S;
  vector[data_agrps] comp_I;
  vector[data_agrps] comp_Im;
  vector[data_agrps] comp_pI;
  
  theta[1] = lambda0;
  ///could just make initial values = y (if works in model running)
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
    //(1) only relevant age groups selected; (2) states calculated only @tmax
    comp_S  = (to_vector(y[t,(data_rows[1:16])])) * tot_pop;  //total no. (age-stratified) in S
    
    comp_I  = (to_vector(y[t,(data_rows[17:32])])) * tot_pop;  //total no. (age-stratified) in I
    
    comp_Im = (to_vector(y[t,(data_rows[33:48])])) * tot_pop;   //total no. (age-stratified) in Im
    
    //compute seroprevalence at tmax
    for(j in 1:data_agrps){  //no. age groups in data
    comp_pI[j] = (comp_I[j] + comp_Im[j]) / (comp_I[j] + comp_Im[j] + comp_S[j]);
    }
}

model {
  lambda0 ~ lognormal(log(.1), .01);
  
  //debug
  if(doprint==1) {
    print("lambda0: ", lambda0);
    //print("S: ", y[t,1:400]);
    print("comp_S: ", comp_S);
    print("comp_I: ", comp_I);
    print("comp_Im: ", comp_Im);
    print("comp_pI: ", comp_pI);
  }
  
  // likelihood 
  if (inference==1) {
    for(j in 1:data_agrps) {  //no. age groups in data
      target += binomial_lpmf(cases[j] | n[j], comp_pI[j]);  //evaluate likelihood at tmax
    }
  }
}
