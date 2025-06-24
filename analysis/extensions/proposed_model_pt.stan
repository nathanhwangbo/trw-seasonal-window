//
// This Stan program defines a simple linear regression model,
// where we model trw as a linear function of observed climate (clim)
//

// The input data is two vectors of length 'N': trw and clim
data {
  int<lower=0> N_years; // total number of years we have ANY tree ring data
  int<lower=0> N_years_obs_precip; // total number of years we have clim data
  int<lower=0> N_years_obs_temp; // total number of years we have clim data
  int<lower=0> N_trees;

  vector[N_years_obs_precip * 12] precip_obs;
  vector[N_years_obs_temp * 12] temp_obs;


  array[N_trees] int s; // length of each tree ring series
  array[N_trees] int start_years; // start year of each tree
  vector[sum(s)] trw; // long vector with all trees concatenated
  vector[sum(s)] age; // long vector with all trees concatenated

  vector[12] w_param_precip; // initial guess of what the weights should look like
  vector[12] w_param_temp; // initial guess of what the weights should look like

  // stuff for generated quantities (not part of the model)
  simplex[12] w_djf;
  vector[12] month_means_mod_p;
  vector[12] month_sds_mod_p;

  simplex[12] w_jja;
  vector[12] month_means_mod_t;
  vector[12] month_sds_mod_t;

}

transformed data {

  int<lower=0> N_months;
  int<lower=0> N_months_mis_precip;
  int<lower=0> N_months_obs_precip;
  int<lower=0> N_years_mis_precip;
  int<lower=0> reconstruction_start_year;
  array[N_trees] int start_index; // the clim index corresponding to the tree's start year

  int<lower=0> N_months_mis_temp;
  int<lower=0> N_months_obs_temp;
  int<lower=0> N_years_mis_temp;


  N_months = N_years * 12;

  N_years_mis_precip = N_years - N_years_obs_precip;
  N_months_mis_precip = N_years_mis_precip * 12;
  N_months_obs_precip = N_years_obs_precip * 12;

  N_years_mis_temp = N_years - N_years_obs_temp;
  N_months_mis_temp = N_years_mis_temp * 12;
  N_months_obs_temp = N_years_obs_temp * 12;

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
  real beta_t;
  real beta_p;
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
  real mu_precip;
  real<lower=0, upper=1> alpha_ar_precip;
  real<lower=0> sigma_ar_precip;

  // autoregressive parameters on monthly climate
  real mu_temp;
  real<lower=0, upper=1> alpha_ar_temp;
  real<lower=0> sigma_ar_temp;


  // the reconstruction
  vector[N_months_mis_precip] precip_mis;
  vector[N_months_mis_temp] temp_mis;


  // parameters for the inferred monthly weights
  // weights: https://statmodeling.stat.columbia.edu/2009/04/29/conjugate_prior/
  vector[12] theta_precip;
  real<lower=0> sq_rho_w_precip;
  real<lower=0> sq_alpha_w_precip;
  vector<lower=0>[12] sq_sigma_w_precip;

  // temp weights
  vector[12] theta_temp;
  real<lower=0> sq_rho_w_temp;
  real<lower=0> sq_alpha_w_temp;
  vector<lower=0>[12] sq_sigma_w_temp;





  // parameters for model discrepancy term
  vector[N_years] eta_eps;
  real<lower=0> sigma_eta;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> alpha2;


}


transformed parameters {
  // concat the missing and observed clim
  vector[N_months] precip_monthly;
  vector[N_years] precip_w_yearly; // weighted monthly clim
  precip_monthly[1:N_months_mis_precip] = precip_mis;
  precip_monthly[(N_months_mis_precip + 1):N_months] = precip_obs;

  vector[N_months] temp_monthly;
  vector[N_years] temp_w_yearly; // weighted monthly clim
  temp_monthly[1:N_months_mis_temp] = temp_mis;
  temp_monthly[(N_months_mis_temp + 1):N_months] = temp_obs;


  // // defining weights
  vector[12] w_p;
  for(i in 1:12){
    w_p[i] = exp(theta_precip[i]) / sum(exp(theta_precip));
  }

  // temperature
  vector[12] w_t;
  for(i in 1:12){
    w_t[i] = exp(theta_temp[i]) / sum(exp(theta_temp));
  }

  { // creating a block -- pos won't be a parameter
    int pos;
    pos = 1;
    for(t in 1:N_years){
      // what we're modelling
      precip_w_yearly[t] = dot_product(w_p, segment(precip_monthly, pos, 12));
      temp_w_yearly[t] = dot_product(w_t, segment(temp_monthly, pos, 12));

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
  beta_t ~ normal(0, 5);
  beta_p ~ normal(0, 5);

  mu_s ~ normal(0, 5);
  logsigma_s ~ normal(0, 1);
  sigma ~ normal(mu_s, sigma_s);

  beta_age ~ normal(0,1);

  mu_precip ~ normal(0,2);
  sigma_ar_precip ~ normal(0, 5);
  alpha_ar_precip ~ std_normal();

  mu_precip ~ normal(0,2);
  sigma_ar_temp ~ normal(0, 5);
  alpha_ar_temp ~ std_normal();

  precip_mis[1] ~ normal(mu_precip,sigma_ar_precip);
  temp_mis[1] ~ normal(mu_temp,sigma_ar_temp);



  // setting up weights covariance

  matrix[12, 12] L_K_precip;
  matrix[12, 12] K_precip;

  sq_rho_w_precip ~ inv_gamma(5,5);
  sq_alpha_w_precip ~ std_normal();
  sq_sigma_w_precip ~ normal(0, 1);

  for (i in 1:(12 - 1)) {
    K_precip[i, i] = sq_alpha_w_precip + sq_sigma_w_precip[i] + 0.001;
    for (j in (i + 1):12) {
      K_precip[i, j] = sq_alpha_w_precip * exp(-1/(2*sq_rho_w_precip) * square(min(j - i, i + 12 - j)));
      K_precip[j, i] = K_precip[i, j];
    }
  }
  K_precip[12, 12] = sq_alpha_w_precip + sq_sigma_w_precip[12] + 0.001;
  L_K_precip = cholesky_decompose(K_precip);
  theta_precip ~ multi_normal_cholesky(w_param_precip, L_K_precip);


  // repeat for temp
  // setting up weights covariance

  matrix[12, 12] L_K_temp;
  matrix[12, 12] K_temp;

  sq_rho_w_temp ~ inv_gamma(5,5);
  sq_alpha_w_temp ~ std_normal();
  sq_sigma_w_temp ~ normal(0, 1);

  for (i in 1:(12 - 1)) {
    K_temp[i, i] = sq_alpha_w_temp + sq_sigma_w_temp[i] + 0.001;
    for (j in (i + 1):12) {
      K_temp[i, j] = sq_alpha_w_temp * exp(-1/(2*sq_rho_w_temp) * square(min(j - i, i + 12 - j)));
      K_temp[j, i] = K_temp[i, j];
    }
  }
  K_temp[12, 12] = sq_alpha_w_temp + sq_sigma_w_temp[12] + 0.001;
  L_K_temp = cholesky_decompose(K_temp);
  theta_temp ~ multi_normal_cholesky(w_param_temp, L_K_temp);



  // likelihood
  int pos;
  pos = 1;
  for(i in 1:N_trees){
    segment(trw, pos, s[i]) ~ normal(beta0[i] +
      beta_p*segment(precip_w_yearly, start_index[i], s[i]) +
      beta_t*segment(temp_w_yearly, start_index[i], s[i]) +
      beta_age[i] * segment(age, pos, s[i]) +
      segment(eta, start_index[i], s[i]), sigma[i]);
    pos = pos + s[i];
  }

  precip_monthly[2:N_months] ~ normal(mu_precip + alpha_ar_precip*(precip_monthly[1:N_months-1] - mu_precip), sigma_ar_precip);
  temp_monthly[2:N_months] ~ normal(mu_temp + alpha_ar_temp*(temp_monthly[1:N_months-1] - mu_temp), sigma_ar_temp);


  alpha0 ~ std_normal();
  alpha1 ~ std_normal();
  alpha2 ~ std_normal();
  sigma_eta ~ normal(0, 0.5);


  eta_eps ~ std_normal();


}


generated quantities {

  vector[N_months] precip_monthly_raw;
  vector[N_years] precip_djf;
  vector[N_months] precip_monthly_log;
  vector[N_years] precip_logdjf;

  vector[N_months] temp_monthly_raw;
  vector[N_years] temp_jja;


  { // creating a block -- pos won't be a parameter
    int pos;
    pos = 1;
    for(t in 1:N_years){

      // get back into un-transformed space
      for(m in 1:12){
        precip_monthly_log[pos + m - 1] = precip_monthly[pos + m - 1]*month_sds_mod_p[m] + month_means_mod_p[m];
        precip_monthly_raw = exp(precip_monthly_log);

        temp_monthly_raw[pos + m - 1] = temp_monthly[pos + m - 1]*month_sds_mod_t[m] + month_means_mod_t[m];


      }

      // what we want: un-transform
      precip_djf[t] = dot_product(w_djf, segment(precip_monthly_raw, pos, 12));
      precip_logdjf[t] = dot_product(w_djf, segment(precip_monthly_log, pos, 12));

      temp_jja[t] = dot_product(w_jja, segment(temp_monthly_raw, pos, 12));



      pos = pos + 12;
    }

  }

}


