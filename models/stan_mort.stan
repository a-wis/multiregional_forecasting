data {
  //LC for mortality
  //
  //Age-sex profile changing
  int<lower=0> R;         // Regions
  int<lower=0> N;         // Age
  int<lower=0> F;         // Forecast horizon
  int<lower=0> T;         // Time
  int<lower=0> S;         // Sex
  matrix [N-1, N-1] mat_R;
  int d[R*N*S*T];   // mortality in vector form
  real p1[R*N*S*T];  // Population
}

transformed data{
  real lp1[R*N*S*T];
  lp1 = log(p1);
}

parameters {
  matrix[N,S] A_a; 
  real G[S];
  matrix[N-1,S] ben;
  real RA[R,N];
  real RS[R,S];
  // real Reg[R];
  real<lower=0> sig1;
  real<lower=0> sig2[3];
  real<lower=0> sig3[2];
  vector<lower=0>[S] sigk;
  
  vector[S] xi[T-1]; //time series
  // real<lower=0, upper=1> slope; 
  vector[S] intercept;
  cholesky_factor_corr[S] L_corr; // prior correlation
  real snd1[R*N*S*T];
  // real snd2[N];
}

transformed parameters {
  matrix[S,S] L_sigma;
  matrix[N,S] A_b;
  real mmd[R*N*S*T]; 
  real mpd[R*N*S*T];
  vector[S] k[T-1]; //time series 

  for (j in 1:S){
    for (i in 1:(N-1)){  
      A_b[i,j] = ben[i,j];
    }
    A_b[N,j] = 1-sum(ben[1:(N-1),j]);
  }
  L_sigma = diag_pre_multiply(sigk, L_corr); 

  for (t in 1:(T-1)) {
    for (s in 1:S) {
      k[t,s] = t*intercept[s] + sum(xi[1:t,s]);
    }
  }

  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:T){
        for (i in 1:R){
          mmd[R*(T*(S*(a-1)+s-1)+t-1)+i] = G[s] + RA[i,a]+RS[i,s] + A_a[a,s] + (t == 1 ? 0 : A_b[a,s]*k[t-1,s]);
        }
      }
    }
  }
  
  for (i in 1:(R*N*T*S)){
    mpd[i] = mmd[i] + sig1*snd1[i] + lp1[i];
  }
  
}


model {
  vector[N-1] mb;
  matrix [N-1, N-1] mat_R1;
  real taub;
  // vector[S] mus[T-1];
  vector[2] vzero = [0.0, 0.0]';
  // priors for precision
  taub = N*N;
  sig1  ~ cauchy(0, 2); 
  sig2  ~ cauchy(0, 2);
  sig3  ~ normal(0, 1);
  
  
  // priors for parameters
  // for (a in 1:N){
    // A[a] ~ normal(0,sig2[1]);
    // }

  for (a in 1:N){
    for (s in 1:S){
      A_a[a,s] ~ normal(0,sig2[s]);
    }
  }
  

  // for (s in 1:S){
    //   G[s] ~ normal(0,sig2[2]);
    // }  
  G ~ normal(0,sig2[3]);
  // Reg ~ normal(0,sig2[3]); //sig3[3]
  // priors for random effects  

  for (i in 1:R){
    for (a in 1:N){
      RA[i,a] ~ normal(0,sig3[1]); //sig3[1]
    }
    for (s in 1:S){
      RS[i,s] ~ normal(0,sig3[2]); //sig3[1]
    }
  }

  // prior for A_b  
  for (i in 1:(N-1)) mb[i] = 1.0/N;
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      mat_R1[i,j] = mat_R[i,j]*taub;
    }
  }
  ben[1:(N-1),1] ~ multi_normal_prec(mb,mat_R1);
  ben[1:(N-1),2] ~ multi_normal_prec(mb,mat_R1);
  


  // time series model AR1
  intercept ~ normal(0,10);
  // slope ~ normal(0,1);
  L_corr ~ lkj_corr_cholesky(2.0);
  sigk ~ normal(0,1);//cauchy(0,2);
  xi ~ multi_normal_cholesky(vzero,L_sigma);
  //Model
  
  snd1 ~ normal(0,1);
  d ~ poisson_log(mpd);
}

generated quantities {
  vector[S] kf[F];
  vector[S] xif[F];
  real mmd_f[R*N*S*F];
  matrix[S,S] Corrm;
  matrix[S,S] Sigma;
  vector[2] vzero = [0.0, 0.0]';
  vector[S] sumk = [sum(xi[1:(T-1),1]), sum(xi[1:(T-1),2])]';
  Sigma = L_sigma * L_sigma';
  Corrm = L_corr * L_corr';

  for (t in 1:F) {
    // musf[t] = intercept + kf[t-1];
    xif[t,1:S] = multi_normal_cholesky_rng(vzero,L_sigma);
  }
  for (s in 1:S) {
    kf[1,s] = T*intercept[s] + xif[1,s] + sumk[s];
  }
  for (t in 2:F) {
    for (s in 1:S) {
      kf[t,s] = (T-1+t)*intercept[s] + sum(xif[1:t,s]) + sumk[s];
    }
  }
  
  // for (i in 1:S) musf[1,i] = intercept[i] + k[T-1,i];  
  // xif[1,1:S] = multi_normal_cholesky_rng(vzero,L_sigma);

  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:F){
        for (i in 1:R){
          mmd_f[R*(F*(S*(a-1)+s-1)+t-1)+i] = G[s] + RA[i,a]+RS[i,s] + A_a[a,s] +  A_b[a,s]*kf[t,s];
        }
      }
    }
  }
}
