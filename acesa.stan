data {
  int<lower=0> T; // number of days
  int<lower=0> P; // number of constituents
  int<lower=0> L; // number of source
  int<lower=0> K; // number of free profile parameters
  matrix[P,T] Xt; //transposed data matrix
}

//hold
parameters {
  vector[K] Fvec;
  matrix<lower=0>[L,T] G; // source contributions
  // matrix<lower=0>[P,L] F; // source profiles
  vector<lower=0>[P] phi2; // variances of X
  vector[L] mu; // means of G
  vector<lower=0>[L] sigma2; // variances of G
}

transformed parameters {
  matrix[P,T] theta; // Mean for x
  theta <- G;
  // matrix<lower=0>[L,L] sigma; //Variance for G
  // sigma <- diag_matrix(sigma2); 
}

model {
  matrix[P,L] F;
  F[row_mark,col_mark] <- Fvec;
  //hold 
  mu ~ normal(0, 1000); // Nikolov et al. 2007
  sigma2 ~ inv_gamma(0.01, 0.01); // Nikolov et al. 2007
  for (l in 1:L)
    G[l] ~ lognormal(mu[l],sigma2[l]);
  phi2 ~ inv_gamma(0.01, 0.01); // Nikolov et al. 2007
  for (p in 1 : P) 
    F[p] ~ lognormal(0, 100000); // ??
  for (p in 1 : P)
    Xt[p] ~ normal(theta[p], phi2[p]);
}


