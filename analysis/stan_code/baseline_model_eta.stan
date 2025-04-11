//
// This Stan program defines a simple linear regression model,
// where we model trw as a linear function of observed climate (clim_yearly)
//
// the goal is to mimic barcast in the time series only setting, and to
// combine it with the hierarchical regression of schofield et al 2015.

// The input data is two vectors of length 'N': trw and clim_yearly
data {
  int<lower=0> N_years; // total number of years we have ANY tree ring data
  int<lower=0> N_years_obs; // total number of years we have clim_yearly data
  int<lower=0> N_trees;

  vector[N_years_obs] clim_obs;

  array[N_trees] int s; // length of each tree ring series
  array[N_trees] int start_years; // start year of each tree
  vector[sum(s)] trw; // long vector with all trees concatenated
  vector[sum(s)] age; // long vector with all ages concatenated

  // for the generated quantities block
  real mean_logdjf;
  real<lower=0> sd_logdjf;
}

transformed data {
  int<lower=0> N_years_mis;
  int<lower=0> reconstruction_start_year;
  array[N_trees] int start_index; // the clim_yearly index corresponding to the tree's start year


  N_years_mis = N_years - N_years_obs;
  reconstruction_start_year = min(start_years);

  for(i in 1:N_trees){
    start_index[i] = start_years[i] - reconstruction_start_year + 1;
  }

}

parameters {
  vector[N_trees] beta0;
  real<lower=0> beta1;
  vector<upper=0>[N_trees] beta_age;
  vector<lower=0>[N_trees] sigma;

  // real<lower=0> sigma_clim;


  real<lower=0> sigma_t;
  real<lower=0, upper=1> alpha;
  real mu;

  vector[N_years_mis] clim_yearly_mis;


  vector[N_years] eta_eps;
  real<lower=0> sigma_eta;


}

transformed parameters {
  vector[N_years] clim_yearly;

  clim_yearly[1:N_years_mis] = clim_yearly_mis;
  clim_yearly[N_years_mis+1:N_years] = clim_obs;


  vector[N_years] eta;
  eta = eta_eps * sigma_eta;


}

model {

  mu ~ normal(0, 2);
  beta0 ~ normal(0, 2);
  beta_age ~ normal(0, 2);
  beta1 ~ normal(0, 2);

  sigma_t ~ normal(0, 2);
  alpha ~ std_normal();


  clim_yearly[1] ~ normal(0, sigma_t);
  clim_yearly[2:N_years] ~ normal(mu + alpha*(clim_yearly[1:N_years-1] - mu), sigma_t);



  sigma_eta ~ normal(0, 2);
  sigma ~ normal(0,2);
  eta_eps ~ std_normal();

  int pos;
  pos = 1;
  for(i in 1:N_trees){
    segment(trw, pos, s[i]) ~ normal(beta0[i] + beta_age[i] * segment(age, pos, s[i]) + beta1 * segment(clim_yearly, start_index[i], s[i]) + segment(eta, start_index[i], s[i]), sigma[i]);
    pos = pos + s[i];
  }





}

generated quantities {
  // note: be mindful of what we pass in as "clim_yearly"
  //    if it's log -> anom, then taking the exp is kinda useless -- it doesn't get us back to raw precip.
  vector[N_years] clim_raw;
  clim_raw = exp(clim_yearly * sd_logdjf + mean_logdjf) - 0.001;
//
//
}
