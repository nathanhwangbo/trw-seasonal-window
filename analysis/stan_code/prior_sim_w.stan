// sampling from the prior for the monthly weights
data {

  vector[12] w_param; // initial guess of what the weights should look like

}

parameters {
  real<lower=0> sq_rho_w;
  real<lower=0> sq_alpha_w;
  vector<lower=0>[12] sq_sigma_w;
}


model {

  sq_rho_w ~ inv_gamma(5, 5);
  sq_alpha_w ~ std_normal();
  sq_sigma_w ~ normal(0, 1);

}

generated quantities {

  matrix[12, 12] K;
  vector[12] theta;
  vector[12] w;

  for (i in 1:(12 - 1)) {
    K[i, i] = sq_alpha_w + sq_sigma_w[i] + 0.001;
    for (j in (i + 1):12) {
      K[i, j] = sq_alpha_w * exp(-1/(2*sq_rho_w) * square(min(j - i, i + 12 - j)));
      K[j, i] = K[i, j];
    }
  }
  K[12, 12] = sq_alpha_w + sq_sigma_w[12] + 0.001;

  theta = multi_normal_rng(w_param, K);

  // normalizing the weights
  for(i in 1:12){
    w[i] = exp(theta[i]) / sum(exp(theta));
  }


}


