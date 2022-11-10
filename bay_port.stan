//https://mc-stan.org/docs/2_23/stan-users-guide/index.html
//https://mc-stan.org/docs/2_23/reference-manual/index.html
data {
  //5 is arbitrary
  int<lower = 5> T;
  int<lower= 2> N;
  real<lower = N - 1> nu;
  real<lower = 0> tau;
  vector[N] eta;
  matrix[T, N] R;
  cov_matrix[N] omega;
}

parameters {
  vector[N] mu;
  cov_matrix[N] sigma;
}

transformed parameters{
  cov_matrix[N] sigma_scaled;
  sigma_scaled = (1 / tau) * sigma;
}

model {
  
  target += inv_wishart_lpdf(sigma | nu, omega);
  target += multi_normal_lpdf(mu | eta, sigma_scaled);
  for(t in 1:T){
    target += multi_normal_lpdf(R[t]| mu, sigma);
  }
  
}

