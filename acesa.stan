data {
  int<lower=0> T; // number of days
  int<lower=0> P; // number of constituents
  int<lower=0> L; // number of source
  int<lower=0> K; // number of free profile parameters
  matrix[P,T] Xt; //transposed data matrix
}

//hold
parameters {
  vector<lower=0>[K] Fvec;
  matrix<lower=0>[L,T] G; // source contributions
  vector<lower=0>[P] phi2; // variances of X
  vector[L] mu; // means of G
  vector<lower=0>[L] sigma2; // variances of G
}


model {

  matrix[P,T] theta; // Mean for x
  matrix[P,L] F; // Source profiles
  for (k in 1:K)
    F[row_mark[k],col_mark[k]] <- Fvec[k];
  //hold 
  mu ~ normal(0, 0.5); // N(0, 1000)Nikolov et al. 2007
  sigma2 ~ inv_gamma(2, 2); // Nikolov et al. 2007 0.01,0.01
  for (l in 1:L)
    G[l] ~ lognormal(mu[l],sigma2[l]);
  phi2 ~ inv_gamma(2, 2); // Nikolov et al. 2007 0.01,0.01
  for (p in 1:P) {
    for (l in 1:L) {  
       F[p,l] ~ normal(0, 100) T[0,]; // ??
      //F[p,l] ~ lognormal(-1, 1) ; // ??
    }
  }
  theta <- F * G;  
  for (p in 1:P)
    Xt[p] ~ normal(theta[p], phi2[p]);
}


