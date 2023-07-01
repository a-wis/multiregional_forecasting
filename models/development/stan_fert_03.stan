data {
  //LC for ifertility
  //
    //Age profile changing, RN, cohort effect
  int<lower=0> R;         // Regions
  int<lower=0> N;         // Age
  int<lower=0> F;         // Forecast horizon
  int<lower=0> T;         // Time
  matrix [N-1, N-1] mat_R;
  matrix [R-1, R-1] mat_RR;
  int d[R*N*T];   // mortality in vector form
  real p1[R*N*T];  // Population
  real sig_t;
}
transformed data{
  real lp1[R*N*T];
  lp1 = log(p1);
}
parameters {
  real cons;
  vector[N] A_a;
  vector[R] R_a;
  vector[N-1] ben;
  vector[R-1] ben2;
  real RA[R,N];
  real<lower=0> sig1;
  real<lower=0> sig2[2];
  real<lower=0> sig3;
  real<lower=0> sigk[2];
  // real<lower=0> sigc;
  vector[T-1] k;
  vector[T-1] k2;
  real intercept[2];
  // real intercept2;
  // vector<lower=0, upper=1>[2] slope;
  // vector[N+T-3+F] coh;
  real<lower=0, upper=1> slope[2];
  real snd1[R*N*T];
  // real snd2[N];
}
transformed parameters {
  vector[N] A_b;
  vector[R] R_b;
  real mmd[R*N*T]; 
  real mpd[R*N*T];
  
  for (i in 1:(N-1)){  
    A_b[i] = ben[i];
  }
  A_b[N] = 1-sum(ben[1:(N-1)]);
  for (i in 1:(R-1)){
    R_b[i] = ben2[i];
    }
  R_b[R] = 1-sum(ben2[1:(R-1)]);

  for (a in 1:N){
    for (t in 1:T){
      for (i in 1:R){
        mmd[R*(T*(a-1)+t-1)+i] = cons + RA[i,a] + A_a[a] + R_a[i] + (t == 1 ? 0 : A_b[a]*k[t-1]) + (t == 1 ? 0 : k2[t-1]*R_b[i]);// + ((a > N-2 && t == 1)) || (a == N && t == 2) ? 0 : coh[N+t-2-a]; 
      }
    }
  }
  
  for (i in 1:(R*N*T)){
    mpd[i] = mmd[i] + sig1*snd1[i] + lp1[i];
  }
  
}
model {
  vector[N-1] mb;
  vector[R-1] mb2;
  matrix[N-1, N-1] mat_R1;
  matrix[R-1, R-1] mat_R2;
  real taub;
  real taub1;
  // priors for precision
  taub = N*N;
  taub1 = R*R;
  
  sig1  ~ student_t(2.5,0,sig_t); 
  // sigk  ~ cauchy(0, 1);
  
  sig2  ~ student_t(2.5,0,sig_t);
  
  sig3  ~ normal(0, 0.2);
  
  
  // priors for parameters
  cons ~ normal(0,5);
  A_a ~ normal(0,sig2[1]);
  R_a ~ normal(0,sig2[2]);

  // priors for random effects  
  for (i in 1:R){
    for (a in 1:N){
      RA[i,a] ~ normal(0,sig3); //sig3[1]
    }
  }
  
  // prior for A_b  
  for (i in 1:(N-1)) mb[i] = 1.0/N;
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      mat_R1[i,j] = mat_R[i,j]*taub;
    }
  }
  for (i in 1:(R-1)) mb2[i] = 1.0/R;
  for (i in 1:(R-1)){
    for (j in 1:(R-1)){
      mat_R2[i,j] = mat_RR[i,j]*taub1;
    }
  }
  ben[1:(N-1)] ~ multi_normal_prec(mb,mat_R1);
  ben2[1:(R-1)] ~ multi_normal_prec(mb2,mat_R2);
  
  // time series model AR1
  intercept ~ normal(0,2);
  // intercept2 ~ normal(0,5);
  slope ~ normal(0.5,0.2);
  sigk ~ normal(0,.2);
  // sigc ~ normal(0,1);
  
  k[1] ~ normal(intercept[1], sigk[1]);
  k[2:(T-1)] ~ normal(intercept[1] + slope[1]*k[1:(T-2)], sigk[1]);
  k2[1] ~ normal(intercept[2], sigk[2]);
  k2[2:(T-1)] ~ normal(intercept[2] + slope[2] * k2[1:(T-2)], sigk[2]);   
  // coh[1] ~ normal(0.0, sigc);
  // coh[2:(N+T-3+F)] ~ normal(slope[2]*coh[1:(N+T-4+F)], sigc);
  
  //Model
  
  snd1 ~ normal(0,1);
  d ~ poisson_log(mpd);
}
generated quantities {
  vector[F] kf;
  vector[F] k2f;
  real mmd_f[R*N*F];

  kf[1] = normal_rng(intercept[1] +slope[1]*k[T-1],sigk[1]);  
  for (t in 2:F) {
    kf[t] = normal_rng(intercept[1] +slope[1]*kf[t-1], sigk[1]);
  }
  k2f[1] = normal_rng(intercept[2] + slope[2] * k2[T-1],sigk[2]);
  for (t in 2:F) {
    k2f[t] = normal_rng(intercept[2] + slope[2] * k2f[t-1],sigk[2]);
  }
  for (a in 1:N){
    for (t in 1:F){
      for (i in 1:R){
        mmd_f[R*(F*(a-1)+t-1)+i] = cons + RA[i,a] + A_a[a] + R_a[i] + A_b[a]*kf[t]+ k2f[t]*R_b[i];// + coh[N+T+t-2-a];
      }
    }
  }
}
