data {
  //LC for ifertility
  //
    //Age profile changing, RN, cohort effect
  int<lower=0> R;         // Regions
  int<lower=0> N;         // Age
  int<lower=0> F;         // Forecast horizon
  int<lower=0> T;         // Time
  matrix [N-1, N-1] mat_R;
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
  // vector[N] A_a;  
  vector[N-1] ben;
  real RA[R,N];
  real<lower=0> sig1;
  real<lower=0> sig2;
  real<lower=0> sig3;
  real<lower=0> sigk;
  real<lower=0> sigc;
  vector[T-1] k;
  // real intercept;
  // real intercept2;
  // vector<lower=0, upper=1>[2] slope;
  vector[N+T-3+F] coh;
  real<lower=0, upper=1> slope[2];
  real snd1[R*N*T];
  // real snd2[N];
}
transformed parameters {
  vector[N] A_b;
  real mmd[R*N*T]; 
  real mpd[R*N*T];
  
  for (i in 1:(N-1)){  
    A_b[i] = ben[i];
  }
  A_b[N] = 1-sum(ben[1:(N-1)]);

  for (a in 1:N){
    for (t in 1:T){
      for (i in 1:R){
        mmd[R*(T*(a-1)+t-1)+i] = cons + RA[i,a] + (t == 1 ? 0 : A_b[a]*k[t-1]) + ((a > N-2 && t == 1)) || (a == N && t == 2) ? 0 : coh[N+t-2-a]; 
      }
    }
  }
  
  for (i in 1:(R*N*T)){
    mpd[i] = mmd[i] + sig1*snd1[i] + lp1[i];
  }
  
}
model {
  vector[N-1] mb;
  matrix[N-1, N-1] mat_R1;
  real taub;
  
  // priors for precision
  taub = N*N;
  sig1  ~ student_t(2.5,0,sig_t); 
  // sigk  ~ cauchy(0, 1);
  
  // sig2  ~ student_t(2.5,0,sig_t); 
  
  sig3  ~ normal(0, 2);
  
  
  // priors for parameters
  // cons ~ normal(0,5);
  // A_a ~ normal(0,sig2);

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
  ben[1:(N-1)] ~ multi_normal_prec(mb,mat_R1);
  
  // time series model AR1
  // intercept ~ normal(0,5);
  // intercept2 ~ normal(0,5);
  slope ~ normal(0.5,0.2);
  sigk ~ normal(0,0.2);
  sigc ~ normal(0,0.2);
  
  k[1] ~ normal(0.0, sigk);
  k[2:(T-1)] ~ normal(slope[1]*k[1:(T-2)], sigk);
    
  coh[1] ~ normal(0.0, sigc);
  coh[2:(N+T-3+F)] ~ normal(slope[2]*coh[1:(N+T-4+F)], sigc);
  
  //Model
  
  snd1 ~ normal(0,1);
  d ~ poisson_log(mpd);
}
generated quantities {
  vector[F] kf;
  real mmd_f[R*N*F];

  kf[1] = normal_rng(slope[1]*k[T-1],sigk);  
  for (t in 2:F) {
    kf[t] = normal_rng(slope[1]*kf[t-1], sigk);
  }
  for (a in 1:N){
    for (t in 1:F){
      for (i in 1:R){
        mmd_f[R*(F*(a-1)+t-1)+i] = cons + RA[i,a] +  A_b[a]*kf[t] + coh[N+T+t-2-a];
      }
    }
  }
}
