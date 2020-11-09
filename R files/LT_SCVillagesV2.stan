//// MCMC :  3 SC VILLAGES - LONG TERM - V2
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
    //Function for DBR close to 0
    real smallDBR (real DBRi){
      if(DBRi>0){
        return DBRi;
      }
      return 0.000001;
    }
    //Function for ODE
    real[] patilleode(real t, real[] y, real[] theta,real[] x_r, int[] x_i){
      //stages
      real E = y[1];   
      real L1 = y[2];   
      real L2 = y[3];
      real L3 = y[4];
      real L4 = y[5];
      real L5 = y[6];
      real L6 = y[7];
      real L7 = y[8];
      real P = y[9];
      real N = y[10];
      real Pa = y[11];
      
      //parameters
      real DBRstar = theta[1]; //from transformed parameters
      real efficacySC = theta[2];
      real decaySC = theta[3];
      //real efficacySC = theta[1];
      //real decaySC = theta[2];
      
      real muEo = x_r[1]; //from transformed data
      real muL = x_r[2];
      real muP = x_r[3];
      real g = x_r[4];
      real muV = x_r[5];
      real DeltaE = x_r[6];
      real DeltaL = x_r[7];
      real DeltaP = x_r[8];
      real HdivHBI = x_r[9];
      real EpsN = x_r[10];
      real EpsP = x_r[11];
      real coef1 = x_r[12];
      real coef2 = x_r[13];
      real coef1bis = x_r[14];
      real coef2bis = x_r[15];
      //real DBRstar = x_r[16];
      
      //derivatives
      real dE_dt;
      real dL1_dt;
      real dL2_dt;
      real dL3_dt;
      real dL4_dt;
      real dL5_dt;
      real dL6_dt;
      real dL7_dt;
      real dP_dt;
      real dN_dt;
      real dPa_dt;
      //local variables
      real mu_sc;
      real mu_rain;
      real betaN;
      real betaP;
      real Vstar;
      real w;
      real K;
      betaN = EpsN*exp(-muV*g)/g;
      betaP = (EpsP * muV) / (exp(muV * g) - 1);
      Vstar = DBRstar * HdivHBI * g;
      w = 2 * DeltaE * (1 + (muL * DeltaL) / 7) ^ 7 * (1 + muP * DeltaP) * (1 + muV *g) * muV;
      K = (Vstar * muEo * (w ^ 2)) / ((1 + g * muV) * (betaN * g * muV + betaP - w /DeltaE - muEo * w));
      //Mortality factors
      if(t>1511){
        mu_sc = -(log(1-efficacySC)/2)*exp(-decaySC*(t-1511));
      }else{
        mu_sc = 0;
      }
      if((t>1536)&&(t<=1612)){
        mu_rain = coef1*sin(2*pi()*(t-1500)/365+0.24)*exp(-coef2*(t-1536));
      }else{
        if((t>1674)&&(t<=1793)){
          mu_rain = coef1bis*-sin(2*pi()*(t-1500)/365+0.24)*exp(-coef2bis*(t-1674));
        }else{
          mu_rain = 0;
        }
      }
      //print("t = ",t);
      //print("mu_rain = ",mu_rain);
      //print("mu_sc = ",mu_sc);
      //ODE
      //rate of change eggs dE/dt
      dE_dt = betaN * N + betaP * Pa   - E * (1 / DeltaE) - (muEo * (1 + E / K) * E) - mu_rain*E - mu_sc*E;
      //rate of change of larval instars 1-7
      dL1_dt = (E * (1 / DeltaE)) - ((7 / DeltaL + muL) * L1) - mu_rain*L1 - mu_sc*L1;
      dL2_dt = 7 * (L1 * (1 / DeltaL)) - ((7 / DeltaL + muL) * L2) - mu_rain*L2 - mu_sc*L2;
      dL3_dt = 7 * (L2 * (1 / DeltaL)) - ((7 / DeltaL + muL) * L3) - mu_rain*L3 - mu_sc*L3;
      dL4_dt = 7 * (L3 * (1 / DeltaL)) - ((7 / DeltaL + muL) * L4) - mu_rain*L4 - mu_sc*L4;
      dL5_dt = 7 * (L4 * (1 / DeltaL)) - ((7 / DeltaL + muL) * L5) - mu_rain*L5 - mu_sc*L5;
      dL6_dt = 7 * (L5 * (1 / DeltaL)) - ((7 / DeltaL + muL) * L6) - mu_rain*L6 - mu_sc*L6;
      dL7_dt = 7 * (L6 * (1 / DeltaL)) - ((7 / DeltaL + muL) * L7) - mu_rain*L7 - mu_sc*L7;
      //rate of change of pupae
      dP_dt = ((7 / DeltaL) * L7) - ((1 / DeltaP + muP) * P) - mu_rain*P - mu_sc*P;
      //rate of change nulliparous adults (female) dN/dt
      dN_dt = 0.5 * P * 1 / DeltaP - (1 / g + muV) * N;
      //rate of change parous adults dPa/dt
      dPa_dt = (1 / g) * N - ((muV) * Pa);
      //return the rate of change             
      return {dE_dt, dL1_dt,dL2_dt, dL3_dt, dL4_dt, dL5_dt, dL6_dt, dL7_dt,dP_dt, dN_dt, dPa_dt};
    }
}
// The input data
data {
  int<lower=1> n_days;  //number of days on which DBR is measured
  real initialpop[11];  //initial population
  real t0;              //t0
  real ts[n_days];      //sample of days
  int DBRs_patille[n_days];     //sample of DBR per day (patille)
  int DBRs_pwomunu[n_days];     //sample of DBR per day (pwomunu)
  int DBRs_eluguB[n_days];     //sample of DBR per day (Elugu A1)
  real muEo;            //constant
  real EpsN;
  real EpsP;
  real <lower=0,upper=0.5>Em;
  real <lower=0,upper=1>muL;
  real <lower=0,upper=1>muP;
  real temp;
  real g;
  real HdivHBI;
  real rainfallmortality1_bajere;
  real rainfalldecay1_bajere;
  real rainfallmortality2_bajere;
  real rainfalldecay2_bajere;
  real rainfallmortality1_okidi;
  real rainfalldecay1_okidi;
  real rainfallmortality2_okidi;
  real rainfalldecay2_okidi;
  real rainfallmortality1_elugua1;
  real rainfalldecay1_elugua1;
  real rainfallmortality2_elugua1;
  real rainfalldecay2_elugua1;
  //real DBRstar_bajere;
  //real DBRstar_okidi;
  //real DBRstar_elugua1;
}
transformed data {
  real water_temp = 0.9844*temp - 1.0352;
  real deltaE = 11.493 * exp(-0.0701 * water_temp);
  real deltaL = 87.527 * exp(-0.0785 * water_temp);
  real deltaP = 20.098 * exp(-0.0699 * water_temp);
  real muV = 0.0027*temp*temp - 0.163*temp + 2.602 + Em;
  //For ODE
  //real x_r_patille[15] = {muEo,muL,muP,g,muV,deltaE,deltaL,deltaP,HdivHBI,EpsN,EpsP,rainfallmortality1_bajere,rainfallmortality2_bajere,rainfalldecay2_bajere,DBRstar_bajere};
  //real x_r_pwomunu[15] = {muEo,muL,muP,g,muV,deltaE,deltaL,deltaP,HdivHBI,EpsN,EpsP,rainfallmortality1_okidi,rainfallmortality2_okidi,rainfalldecay2_okidi,DBRstar_okidi};
  //real x_r_eluguB[15] = {muEo,muL,muP,g,muV,deltaE,deltaL,deltaP,HdivHBI,EpsN,EpsP,rainfallmortality1_elugua1,rainfallmortality2_elugua1,rainfalldecay2_elugua1, DBRstar_elugua1};
  real x_r_patille[15] = {muEo,muL,muP,g,muV,deltaE,deltaL,deltaP,HdivHBI,EpsN,EpsP,rainfallmortality1_bajere,rainfalldecay1_bajere,rainfallmortality2_bajere,rainfalldecay2_bajere};
  real x_r_pwomunu[15] = {muEo,muL,muP,g,muV,deltaE,deltaL,deltaP,HdivHBI,EpsN,EpsP,rainfallmortality1_okidi,rainfalldecay1_okidi,rainfallmortality2_okidi,rainfalldecay2_okidi};
  real x_r_eluguB[15] = {muEo,muL,muP,g,muV,deltaE,deltaL,deltaP,HdivHBI,EpsN,EpsP,rainfallmortality1_elugua1,rainfalldecay1_elugua1,rainfallmortality2_elugua1,rainfalldecay2_elugua1};

  int x_i[0];
}
// The parameters accepted by the model
parameters {
  //real <lower=0,upper=1>efficacySC;
  //real <lower=0>decaySC;
  real <lower = 0,upper = 0.8>efficacySC;
  real <lower= 0,upper= 1>decaySC;
  real <lower=70,upper=400>DBRstar_pwomunu;
  real <lower=70,upper=400>DBRstar_patille;
  real <lower=70,upper=400>DBRstar_eluguB;
}
transformed parameters {
  vector[n_days] lambda_DBR_patille;
  vector[n_days] lambda_DBR_pwomunu;
  vector[n_days] lambda_DBR_eluguB;
  real solvepop_patille[n_days, 11]; //11 stages
  real solvepop_pwomunu[n_days, 11]; //11 stages
  real solvepop_eluguB[n_days, 11]; //11 stages
  print("Log density TF1 ",target()); //(if only one core is used during the fit, the log density is printed)
  {
    real theta_patille[3] = {DBRstar_patille,efficacySC,decaySC};
    real theta_pwomunu[3] = {DBRstar_pwomunu,efficacySC,decaySC};
    real theta_eluguB[3] = {DBRstar_eluguB,efficacySC,decaySC};
    //real theta[2] = {efficacySC,decaySC};
    
    // see what the model output looks like before taking further [i.e. try to return solvepop_patille]
    solvepop_patille = integrate_ode_rk45(patilleode, initialpop, t0, ts, theta_patille, x_r_patille, x_i);  //return matrix (state variables at diff times)
    solvepop_pwomunu = integrate_ode_rk45(patilleode, initialpop, t0, ts, theta_pwomunu, x_r_pwomunu, x_i);
    solvepop_eluguB = integrate_ode_rk45(patilleode, initialpop, t0, ts, theta_eluguB, x_r_eluguB, x_i);
    print("solve!");
  }
  lambda_DBR_patille = (col(to_matrix(solvepop_patille), 10)+col(to_matrix(solvepop_patille), 11)) * (1 / HdivHBI) * (1 / g); //extract particular part of solved model matrix [that want to fit to data]
  lambda_DBR_pwomunu = (col(to_matrix(solvepop_pwomunu), 10)+col(to_matrix(solvepop_pwomunu), 11)) * (1 / HdivHBI) * (1 / g);
  lambda_DBR_eluguB = (col(to_matrix(solvepop_eluguB), 10)+col(to_matrix(solvepop_eluguB), 11)) * (1 / HdivHBI) * (1 / g);
  //print("lambda_DBR =",lambda_DBR_patille[22]);
  for(i in 1:n_days){
    lambda_DBR_patille[i] = smallDBR(lambda_DBR_patille[i]);
    lambda_DBR_pwomunu[i] = smallDBR(lambda_DBR_pwomunu[i]);
    lambda_DBR_eluguB[i] = smallDBR(lambda_DBR_eluguB[i]);
  }
  //print("lambda_DBR =",lambda_DBR_patille[22]);
  print("Log density TF2 ",target());
}
// The model to be estimated. 
model {
  //priors
  DBRstar_patille ~ normal(300,30);  //flat prior - 'uniform'
  DBRstar_pwomunu ~ normal(300,30);
  DBRstar_eluguB ~ normal(300,30);
  efficacySC ~ normal(0.6,0.2);
  decaySC ~ normal(0.2, 0.1);;
  print("Log density M1 ",target());
  //sampling distribution
  DBRs_patille ~ poisson(lambda_DBR_patille); //fitting relevant part of model to data
  DBRs_pwomunu ~ poisson(lambda_DBR_pwomunu);
  DBRs_eluguB ~ poisson(lambda_DBR_eluguB);
  print("Log density M2 ",target());
}
//For the posterior predictive checks
generated quantities { 
  int pred_DBR_patille[n_days];
  int pred_DBR_pwomunu[n_days];
  int pred_DBR_eluguB[n_days];
  pred_DBR_patille = poisson_rng(lambda_DBR_patille);
  pred_DBR_pwomunu = poisson_rng(lambda_DBR_pwomunu);
  pred_DBR_eluguB = poisson_rng(lambda_DBR_eluguB);
}







