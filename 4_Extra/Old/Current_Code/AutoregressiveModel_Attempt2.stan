// This script is based off of "Conor_StanForum_Model.stan", which is a 1D version of this script (which works, unlike this script)

functions{
  real[] m_dts(int N, real[] Z, real[] theta){
    real P[N];
    real H[N];
    P[1] = Z[1]; // Where Z are estimated abundances
    H[1] = Z[2];

    for (t in 1:(N-1)){
      P[t] = y[t-1,1]*theta[1] - y[t-1,2]*(theta[2]*y[t-1,1]/(1 + theta[2]*theta[3]*y[t-1,1])); // where y[t-1,1] means the t-1 position in the first column of data (in this case, P); syntactically, there's something wrong here
      H[t] = theta[4] + y[t-1,2]*(theta[5]*(theta[2]*y[t-1,1]/(1 + theta[2]*theta[3]*y[t-1,1]))-theta[6]);
    }
    return P[];
    return H[];
  }
}

data {
  int<lower=0>N; 
  real<lower=0>y[N,2]; // Where y are population abundances from the data, which has 2 columns
}

//transformed data{
//}

parameters {
  real<lower=0>r; // All parms lower-bounded at 0; O and u upper-bounded at 1
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
  real<lower=0>y_init_1[1]; // Initial values estimated separately
  real<lower=0>y_init_2[1];
  real<lower=0>sigma[2]; // Error
}

transformed parameters{
  real Z[N, 2];
  real theta[6];

  theta[1] = r;
  theta[2] = O;
  theta[3] = h;
  theta[4] = b;
  theta[5] = c;
  theta[6] = u;
  
  Z[1] = m_dts(N,y_init_1,theta);
  Z[2] = m_dts(N,y_init_2,theta);
}

model {
  r~normal(2.5,1); //Priors for parms
  O~normal(0.01,2);
  h~normal(0.07,2);
  b~normal(35,1);
  c~normal(0.3,1);
  u~normal(0.4,1);
  y_init_1~lognormal(log(140), 1); //Priors for inits
  y_init_2~lognormal(log(140), 1);
  sigma~lognormal(-1, 1);
  y[1] ~ lognormal(log(Z[1]), sigma[1]); // Dist of abundance estimates
  y[2] ~ lognormal(log(Z[2]), sigma[2]);
}

generated quantities { // y_rep are the abundance estimates that you can extract from the output
  real y_rep[N, 2];
  y_rep[1] = lognormal_rng(log(Z[1]), sigma[1]);
  y_rep[2] = lognormal_rng(log(Z[2]), sigma[2]);
}
