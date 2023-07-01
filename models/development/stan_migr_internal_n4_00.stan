data {
  //building upon LC model 
  //Poisson 
  //Age profile changing, OD profile changing
  //S+OA+DA
  //RW
  int<lower=0> R;         // origins and destinations (R x R-1)
  int<lower=0> Co;        // Corridors (R x R-1)
  int<lower=0> N;         // Age
  int<lower=0> F;         // Forecast horizon
  int<lower=0> T;         // Time
  int<lower=0> S;         // Sex
  matrix [Co-1, Co-1] mat_Rc;
  matrix [N-1, N-1] mat_Ra;
  int d[Co*N*S*T];   // OD migration in vector form
  real p1[Co*N*S*T];  // Population
  int corridor[Co];       // Corridors index for destinations
  int unicorr[Co];       // uni-directional corridors index 
  int ind_DA[Co-R+1];
  int ind_DA0[R-1];
  // prior
  real sig_t;
}

transformed data{
  int Co2;
  real lp1[Co*N*S*T];
  Co2 = Co/2;
  lp1 = log(p1);
}

parameters {
  vector[Co2] OD_uni;
  real G;
  vector[N-1] ben1;
  vector[Co-1] ben2;
  real OA[R-1,N-1];
  real DA[R-1,N-1];
  real<lower=0> sig1;
  real<lower=0> sig2;
  real<lower=0> sig3[4];
  real<lower=0> sigk[2];
  vector[T-1] k1; //time series
  vector[T-1] k2; //time series 
  // real<lower=0, upper=1> slope[2]; 
  // real intercept[2];
  real snd1[Co*N*S*T];
  // real Grnd;
  real A_arnd[N];
  real odau_snd[Co];
  // real snd2[N];
}

transformed parameters {
  vector[N] A_a;
  vector[N] A_b;
  vector[Co] OD_b;
  real mmd[Co*N*S*T]; 
  real mpd[Co*N*S*T];
  real cOA[Co,N];
  real cDA[Co,N];
  vector[Co] OD_a;

  
  for (i in 1:N) A_a[i] = A_arnd[i] * sig2;
  // G = sig2[2]*Grnd;
  for (r in 1:Co) OD_a[r] = OD_uni[unicorr[r]] + sig3[3]*odau_snd[r];
  
  for (i in 1:(N-1)){  
    A_b[i] = ben1[i];
  }
  A_b[N] = 1-sum(ben1[1:(N-1)]);
  
  for (i in 1:(Co-1)){  
    OD_b[i] = ben2[i];
  }
  OD_b[Co] = 1-sum(ben2[1:(Co-1)]);
  
  for (i in 1:(R-1)){
    for (j in 1:(R-1)){
      cOA[(R-1)*(i-1)+j,1] = 0;
      for (a in 2:N){
        cOA[(R-1)*(i-1)+j,a] = OA[i,a-1]; 
      }
    }
  }
  
  for (j in 1:(R-1)){
    for (a in 1:N){
      cOA[(R-1)*(R-1)+j,a] = 0; 
    }
  }
  
  
  for (i in ind_DA){
    cDA[i,1] = 0;
    for (a in 2:N){
      cDA[i,a] = DA[corridor[i],a-1];
    }
  }
  for (i in ind_DA0){
    for (a in 1:N){
      cDA[i,a] = 0;
    }
  }
  
  
  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:T){
        for (i in 1:Co){
          mmd[Co*(T*(S*(a-1)+s-1)+t-1)+i] = G*(s-1) + OD_a[i] + cOA[i,a] + cDA[i,a]+ A_a[a] + (t == 1 ? 0 : A_b[a]*k1[t-1]) + (t == 1 ? 0 : OD_b[i]*k2[t-1]);
        }
      }
    }
  }
  
  for (i in 1:(Co*N*T*S)){
    mpd[i] = mmd[i] + sig1*snd1[i] + lp1[i];
  }  
  
}


model {
  vector[N-1] mb1;
  vector[Co-1] mb2;
  matrix [N-1, N-1] mat_R1;
  matrix [Co-1, Co-1] mat_R2;
  real taub1;
  real taub2;
  
  // priors for precision
  taub1 = N*N;
  taub2 = Co*Co;
  sig1  ~ student_t(2.5,0,sig_t);
  sigk[1]  ~ student_t(2.5,0,sig_t);
  sigk[2]  ~ student_t(2.5,0,sig_t);
  sig2  ~ student_t(2.5,0,sig_t);
  // for (j in 1:2){  
  //   sig2[j]  ~ student_t(2.5,0,25);
  // }
  for (j in 1:4){  
    sig3[j]  ~ normal(0, 5);
  }
  
  // priors for parameters
  // for (a in 1:N){
    // A[a] ~ normal(0,sig2[1]);
    // }
  
  // A_a ~ normal(0,sig2[1]);
  
  

  
  for (r in 1:Co2){
    OD_uni[r] ~ normal(0,sig3[4]);
  }
  
  // for (s in 1:S){
    //   G[s] ~ normal(0,sig2[2]);
    // }  
  G ~ normal(0,1);
  
  // priors for random effects  
  for (i in 1:(R-1)){
    for (a in 1:(N-1)){
      OA[i,a] ~ normal(0,sig3[1]); //sig3[1]
      DA[i,a] ~ normal(0,sig3[2]); //sig3[2]
    }
  }
  
  
  
  // prior for A_b  
  for (i in 1:(N-1)) mb1[i] = 1.0/N;
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      mat_R1[i,j] = mat_Ra[i,j]*taub1;
    }
  }
  ben1 ~ multi_normal_prec(mb1,mat_R1);
  
  // prior for OD_b  
  for (i in 1:(Co-1)) mb2[i] = 1.0/Co;
  for (i in 1:(Co-1)){
    for (j in 1:(Co-1)){
      mat_R2[i,j] = mat_Rc[i,j]*taub2;
    }
  }
  ben2 ~ multi_normal_prec(mb2,mat_R2);
  
  // time series model AR1
  // intercept ~ normal(0,10);
  // slope ~ normal(0,1);
  
  k1[1] ~ normal( 0.0, sigk[1]);
  k1[2:(T-1)] ~ normal(k1[1:(T-2)], sigk[1]);
  k2[1] ~ normal( 0.0, sigk[2]);
  k2[2:(T-1)] ~ normal(k2[1:(T-2)], sigk[2]);
  
  //Model
  snd1 ~ normal(0,1);
  // Grnd ~ normal(0,1);
  A_arnd ~ normal(0,1);
  odau_snd ~ normal(0,1);
  d ~ poisson_log(mpd);
  
}

generated quantities {
  real kf1[F];
  real kf2[F];
  real mmd_f[Co*N*S*F];
  vector[Co*N*S*T] log_lik;
   
  kf1[1] = normal_rng( k1[T-1], sigk[1]);
  kf2[1] = normal_rng( k2[T-1], sigk[2]);
  for (t in 2:F) kf1[t] = normal_rng(kf1[t-1], sigk[1]);
  for (t in 2:F) kf2[t] = normal_rng(kf2[t-1], sigk[2]); 
  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:F){
        for (i in 1:Co){
          mmd_f[Co*(F*(S*(a-1)+s-1)+t-1)+i] = G*(s-1) + OD_a[i] + cOA[i,a] + cDA[i,a]+A_a[a] + A_b[a]*kf1[t] + OD_b[i]*kf2[t];
        }
      }
    }
  }
  
  for (i in 1:(Co*N*S*T)) log_lik[i] = poisson_log_lpmf(d[i]|mpd[i]);
  
}
