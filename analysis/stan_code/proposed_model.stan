//
// This Stan program defines a simple linear regression model,
// where we model trw as a linear function of observed climate (clim)
//

// The input data is two vectors of length 'N': trw and clim
data {
  int<lower=0> N_years; // total number of years we have ANY tree ring data
  int<lower=0> N_years_obs; // total number of years we have clim data
  int<lower=0> N_trees;

  vector[N_years_obs * 12] clim_obs;
  array[N_trees] int s; // length of each tree ring series
  array[N_trees] int start_years; // start year of each tree
  vector[sum(s)] trw; // long vector with all trees concatenated
  vector[sum(s)] age; // long vector with all trees concatenated

  vector[12] w_param; // initial guess of what the weights should look like

  // stuff for generated quantities (not part of the model)
  simplex[12] w_djf;
  vector[12] month_means_mod;
  vector[12] month_sds_mod;

}

transformed data {

  int<lower=0> N_months;
  int<lower=0> N_months_mis;
  int<lower=0> N_months_obs;

  int<lower=0> N_years_mis;
  int<lower=0> reconstruction_start_year;
  array[N_trees] int start_index; // the clim index corresponding to the tree's start year

  N_months = N_years * 12;
  N_years_mis = N_years - N_years_obs;
  N_months_mis = N_years_mis * 12;
  N_months_obs = N_years_obs * 12;

  // stacking all of the trees into a single vector.
  // we need to keep track of where the splits are.
  reconstruction_start_year = min(start_years);
  for(i in 1:N_trees){
    start_index[i] = start_years[i] - reconstruction_start_year + 1;
  }

  // small helper to keep track of which month is which
  array[12] real month_num;
  for(m in 1:12){
    month_num[m] = m;
  }


}

parameters {

  // calibration regression coefficients
  vector[N_trees] beta0_raw;
  // real<lower=0> beta1;
  real beta1;
  // real beta_age;
  vector[N_trees] beta_age;

  vector<lower=0>[N_trees] sigma; // calibration noise

  // hyper-priors on calibration intercept
  real mu_b0;
  real<lower=0> sigma_b0;

  // hyper-priors on calibration variance (sigma)
  real<lower=0> mu_s;
  real logsigma_s;


  // autoregressive parameters on monthly climate
  real mu;
  real<lower=0, upper=1> alpha_ar;
  real<lower=0> sigma_ar;


  // the reconstruction
  vector[N_months_mis] clim_mis;


  // parameters for the inferred monthly weights
  // weights: https://statmodeling.stat.columbia.edu/2009/04/29/conjugate_prior/
  vector[12] theta;
  real<lower=0> sq_rho_w;
  real<lower=0> sq_alpha_w;
  vector<lower=0>[12] sq_sigma_w;


  // parameters for model discrepancy term
  vector[N_years] eta_eps;
  real<lower=0> sigma_eta;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> alpha2;







}


transformed parameters {
  // concat the missing and observed clim
  vector[N_months] clim_monthly;
  vector[N_years] clim_w_yearly; // weighted monthly clim


  clim_monthly[1:N_months_mis] = clim_mis;
  clim_monthly[(N_months_mis + 1):N_months] = clim_obs;



  // // defining weights
  vector[12] w;
  for(i in 1:12){
    w[i] = exp(theta[i]) / sum(exp(theta));
  }



  { // creating a block -- pos won't be a parameter
    int pos;
    pos = 1;
    for(t in 1:N_years){
      // what we're modelling
      clim_w_yearly[t] = dot_product(w, segment(clim_monthly, pos, 12));



      pos = pos + 12;
    }
  }




  vector[N_years] eta;

  array[N_years] real<lower=0> sigma_eta_vec;
  sigma_eta_vec[1] = sigma_eta;
  eta[1] = eta_eps[1] * sigma_eta;
  for(t in 2:N_years){
    sigma_eta_vec[t] = sqrt(alpha0 + alpha1 * square(eta[t-1]) + alpha2 * square(sigma_eta_vec[t-1]));
    eta[t] = eta_eps[t] * sigma_eta_vec[t];
  }





  vector[N_trees] beta0;
  beta0 = mu_b0 + sigma_b0 * beta0_raw;

  real<lower=0> sigma_s;
  sigma_s = exp(logsigma_s);


}

model {
  mu_b0 ~ normal(0, 5);
  sigma_b0 ~ normal(0,5);
  beta0_raw ~ std_normal();
  // beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);

  mu_s ~ normal(0, 5);
  logsigma_s ~ normal(0, 1);
  sigma ~ normal(mu_s, sigma_s);

  beta_age ~ normal(0,1);

  mu ~ normal(0,2);

  sigma_ar ~ normal(.5, .5);
  alpha_ar ~ std_normal();

  clim_mis[1] ~ normal(mu,sigma_ar);



  // setting up weights covariance

  matrix[12, 12] L_K;
  matrix[12, 12] K;

  sq_rho_w ~ inv_gamma(5,5);
  sq_alpha_w ~ std_normal();
  sq_sigma_w ~ normal(0, 1);

  for (i in 1:(12 - 1)) {
    K[i, i] = sq_alpha_w + sq_sigma_w[i] + 0.001;
    for (j in (i + 1):12) {
      K[i, j] = sq_alpha_w * exp(-1/(2*sq_rho_w) * square(min(j - i, i + 12 - j)));
      K[j, i] = K[i, j];
    }
  }
  K[12, 12] = sq_alpha_w + sq_sigma_w[12] + 0.001;
  L_K = cholesky_decompose(K);
  theta ~ multi_normal_cholesky(w_param, L_K);





  int pos;
  pos = 1;
  for(i in 1:N_trees){
    segment(trw, pos, s[i]) ~ normal(beta0[i] +  beta1*segment(clim_w_yearly, start_index[i], s[i]) + beta_age[i] * segment(age, pos, s[i]) + segment(eta, start_index[i], s[i]), sigma[i]);
    pos = pos + s[i];
  }

  clim_monthly[2:N_months] ~ normal(mu + alpha_ar*(clim_monthly[1:N_months-1] - mu), sigma_ar);


  alpha0 ~ std_normal();
  alpha1 ~ std_normal();
  alpha2 ~ std_normal();
  sigma_eta ~ normal(0, 0.5);


  eta_eps ~ std_normal();


}


generated quantities {

  vector[N_months] clim_monthly_raw;
  vector[N_years] clim_djf;
  vector[N_months] clim_monthly_log;
  vector[N_years] clim_logdjf;

  { // creating a block -- pos won't be a parameter
    int pos;
    pos = 1;
    for(t in 1:N_years){

      // get back into un-transformed space
      for(m in 1:12){
        clim_monthly_log[pos + m - 1] = clim_monthly[pos + m - 1]*month_sds_mod[m] + month_means_mod[m];
        clim_monthly_raw = exp(clim_monthly_log);

      }

      // what we want: un-transform
      clim_djf[t] = dot_product(w_djf, segment(clim_monthly_raw, pos, 12));
      clim_logdjf[t] = dot_product(w_djf, segment(clim_monthly_log, pos, 12));


      pos = pos + 12;
    }

  }

}


