// logisitc dts model
functions{
  real[] m_dts(int T, real[] z, real[] theta){
    
    real n[T];
    
    n[1] = z[1];
    
    for(t in 1:(T-1)){
      n[t+1] = n[t] + n[t] * theta[1] * (1 - n[t]/theta[2]);
    }
    
    return n[];
  }
}

data{
  int<lower=0> T;
  real y[T];
}

transformed data{
}

parameters{
  real<lower=0> y_init[1];
  real<lower=0> r;
  real<lower=0> k;
  real<lower=0> sigma;
}

transformed parameters{
  real z[T];
  real theta[2];
  
  theta[1] = r;
  theta[2] = k  * 100;
  
  z = m_dts(T, y_init, theta);
}

model{
  k ~ normal(0, 1);
  r ~ normal(0, 1);
  y_init ~ normal(0, 1);
  sigma ~ lognormal(-1, 1);
  y ~ lognormal(log(z), sigma);
}

generated quantities{
  real z_pred[T];
  z_pred = lognormal_rng(log(z), sigma);
}
