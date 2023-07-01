data {
  //LC for int immigration
  //
    //Age-sex profile changing, RT + AST, definition effect
    //t priors rather than normal for sig2 in 004; def is random
  int<lower=0> R;         // Regions
  int<lower=0> N;         // Age
  int<lower=0> F;         // Forecast horizon
  int<lower=0> T;         // Time
  int<lower=0> S;         // Sex
  matrix [N-1, N-1] mat_RA;
  // matrix [R-1, R-1] mat_RR;
  int d[R*N*S*T];   // mortality in vector form
  real p1[R*N*S*T];  // Population
  real sig_t;
}
transformed data{
  real lp1[R*N*S*T];
  lp1 = log(p1);
}
parameters {
  real cons;
  matrix[N,S] A_a;
  // real R_a[R];
  // real G;
  matrix[N-1,S] ben;
  // vector[R-1] ben2;
  real RA[R,N];
  real RS[R,S];
  // real Reg[R];
  real<lower=0> sig1;
  real<lower=0> sig2;
  real<lower=0> sig3[2];
  // vector<lower=0>[S] sigk;
  // real<lower=0> sigk2;
  vector[S] k[T-1];
  // vector[T-1] k2;
  vector[S] intercept;
  // real intercept_m;
  cholesky_factor_corr[S] L_corr; // prior correlation
  real snd1[R*N*S*T];
  // real snd2[N];
}
transformed parameters {
  matrix[S,S] L_sigma;
  matrix[N,S] A_b;
  real mmd[R*N*S*T];
  real mpd[R*N*S*T];

  for (j in 1:S){
    for (i in 1:(N-1)){
      A_b[i,j] = ben[i,j];
    }
    A_b[N,j] = 1-sum(ben[1:(N-1),j]);
  }
  L_sigma = diag_pre_multiply([0.1, 0.1]', L_corr);

  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:T){
        for (i in 1:R){
          mmd[R*(T*(S*(a-1)+s-1)+t-1)+i] = cons + RA[i,a] + A_a[a,s] + (t == 1 ? 0 : A_b[a,s]*k[t-1,s]) + RS[i,s];
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
  matrix[N-1, N-1] mat_R1;
  real taub;
  vector[S] mus[T-1];

  // priors for precision
  taub = N*N;
  sig1  ~ student_t(2.5,0,sig_t);
  // sigk  ~ cauchy(0, 1);

  sig2  ~ student_t(2.5,0,sig_t);

  sig3  ~ normal(0, 0.5);
  cons ~ normal(0,5);

  // priors for parameters
  for (a in 1:N){
    for (s in 1:S){
      A_a[a,s] ~ normal(0,sig2);
    }
  }

  // priors for random effects
  // G ~ normal(0,sig2[3]); //sig3[3]
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
      mat_R1[i,j] = mat_RA[i,j]*taub;
    }
  }
  ben[1:(N-1),1] ~ multi_normal_prec(mb,mat_R1);
  ben[1:(N-1),2] ~ multi_normal_prec(mb,mat_R1);
  // time series model AR1
  intercept ~ normal(0,5);
  // intercept_m ~ normal(0,5);
  // slope ~ normal(0.5,0.2);
  L_corr ~ lkj_corr_cholesky(2.0);
  // sigk ~ normal(0,1);

  mus[1] = intercept; //
  for (t in 2:(T-1)) {
    for (i in 1:S) {
      mus[t,i] = intercept[i] + k[t-1,i];
    }
  }
  k ~ multi_normal_cholesky(mus,L_sigma);
  //Model

  snd1 ~ normal(0,1);
  d ~ poisson_log(mpd);
}
generated quantities {
  vector[S] kf[F];
  vector[S] musf[F];
  real mmd_f[R*N*S*F];
  matrix[S,S] Corrm;
  matrix[S,S] Sigma;
  Sigma = L_sigma * L_sigma';
  Corrm = L_corr * L_corr';

  for (i in 1:S) musf[1,i] = intercept[i] + k[T-1,i];
  kf[1,1:S] = multi_normal_cholesky_rng(musf[1,1:S],L_sigma);
  for (t in 2:F) {
    for (i in 1:S) {
      musf[t,i] = intercept[i] + kf[t-1,i];
    }
    kf[t,1:S] = multi_normal_cholesky_rng(musf[t,1:S],L_sigma);
  }
  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:F){
        for (i in 1:R){
          mmd_f[R*(F*(S*(a-1)+s-1)+t-1)+i] = cons + RA[i,a] + A_a[a,s] +  A_b[a,s]*kf[t,s] + RS[i,s];
        }
      }
    }
  }
} 