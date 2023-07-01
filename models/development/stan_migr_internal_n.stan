data {
  //building upon LC model 
  //Poisson 
  //Age profile changing, OD profile changing
  int<lower=0> Co;        // Corridors (R x R-1)
  int<lower=0> N;         // Number of observations: Co*A*S*T
  int<lower=0> F;         // Forecast horizon:Co*A*S*F
  int<lower=0> K;         // number of effects
  int<lower=0> K1;         // number of bilinear effects
  matrix [Co-1, Co-1] mat_R;
  int d[N];   // OD migration in vector form
  real p1[N];  // Population
  matrix[N, K] X; //design matrix
  matrix[N, K1] X1; //design matrix
  int corridor[Co];       // Corridors index for destinations
  int unicorr[Co];       // uni-directional corridors index 
}

transformed data{
  int Co2;
  vector[N] lp1 = log(p1);
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  Co2 = Co/2;
  
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1]; //- means_X[i - 1];
  }
}

parameters {
  vector[N] A_a; 
  vector[Co] OD_a;
  vector[Co2] OD_uni;
  real G;
  vector[N-1] ben1;
  vector[Co-1] ben2;
  real OA[R,N];
  real DA[R,N];
  real O[R-1];
  real D[R];
  real<lower=0> sig1;
  real<lower=0> sig2[4];
  real<lower=0> sig3[4];
  real<lower=0> sigk[2];
  vector[T-1] k1; //time series
  vector[T-1] k2; //time series 
  real<lower=0, upper=1> slope[2]; 
  real intercept[2];
  real snd1[Co*N*S*T];
  // real snd2[N];
}

transformed parameters {
  
  vector[N] A_b;
  vector[Co] OD_b;
  real mmd[Co*N*S*T]; 
  real mpd[Co*N*S*T];
  real cOA[Co,N];
  real cDA[Co,N];
  real cO[Co];
  real cD[Co];
  
  for (i in 1:(N-1)){  
    A_b[i] = ben1[i];
  }
  A_b[N] = 1-sum(ben1[1:(N-1)]);
  
  for (i in 1:(Co-1)){  
    OD_b[i] = ben2[i];
  }
  OD_b[Co] = 1-sum(ben2[1:(Co-1)]);
  
  for (i in 1:R){
    for (j in 1:(R-1)){
      for (a in 1:N){
        cOA[(R-1)*(i-1)+j,a] = OA[i,a]; 
      }
      cO[(R-1)*(i-1)+j] = (i == 1 ? 0 : O[i-1]); 
    }
  }
  
  for (i in 1:Co){
    for (a in 1:N){
      cDA[i,a] = DA[corridor[i],a];
    }
    cD[i] = D[corridor[i]];
  }
  
  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:T){
        for (i in 1:Co){
          mmd[Co*(T*(S*(a-1)+s-1)+t-1)+i] = G*(s-1) + OD_a[i] + cOA[i,a]+cDA[i,a]+cO[i]+cD[i]+ A_a[a] + (t == 1 ? 0 : A_b[a]*k1[t-1]) + (t == 1 ? 0 : OD_b[i]*k2[t-1]);
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
  matrix [N-1, N-1] mat_R1;
  vector[Co-1] mb2;
  matrix [Co-1, Co-1] mat_R2;
  real taub1;
  real taub2;
  
  // priors for precision
  taub1 = N*N;
  taub2 = Co*Co;
  sig1  ~ cauchy(0, 2); 
  sigk[1]  ~ cauchy(0, 2);
  sigk[2]  ~ cauchy(0, 2);
  for (j in 1:4){  
    sig2[j]  ~ cauchy(0, 2);
  }
  for (j in 1:4){  
    sig3[j]  ~ normal(0, 1);
  }
  
  // priors for parameters
  // for (a in 1:N){
    // A[a] ~ normal(0,sig2[1]);
    // }
  
  A_a ~ normal(0,sig2[1]);
  
  
  for (r in 1:Co){
    OD_a[r] ~ normal(OD_uni[unicorr[r]],sig3[3]);
  }
  
  for (r in 1:Co2){
    OD_uni[r] ~ normal(0,sig3[4]);
  }
  
  // for (s in 1:S){
    //   G[s] ~ normal(0,sig2[2]);
    // }  
  G ~ normal(0,sig2[2]);
  
  // priors for random effects  
  O ~ normal(0,sig2[3]); //sig3[3]
  D ~ normal(0,sig2[4]); //sig3[4]
  for (i in 1:R){
    for (a in 1:N){
      OA[i,a] ~ normal(0,sig3[1]); //sig3[1]
      DA[i,a] ~ normal(0,sig3[2]); //sig3[2]
    }
  }
  
  
  
  // prior for A_b  
  for (i in 1:(N-1)) mb1[i] = 1.0/N;
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      mat_R1[i,j] = mat_R[i,j]*taub1;
    }
  }
  ben1 ~ multi_normal_prec(mb1,mat_R1);
  
  // prior for OD_b  
  for (i in 1:(Co-1)) mb2[i] = 1.0/Co;
  for (i in 1:(Co-1)){
    for (j in 1:(Co-1)){
      mat_R2[i,j] = mat_R[i,j]*taub2;
    }
  }
  ben2 ~ multi_normal_prec(mb2,mat_R2);
  
  // time series model AR1
  intercept ~ normal(0,10);
  slope ~ normal(0,1);
  
  k1[1] ~ normal( intercept[1], sigk[1]);
  k1[2:(T-1)] ~ normal(intercept[1] + slope[1]*k1[1:(T-2)], sigk[1]);
  k2[1] ~ normal( intercept[2], sigk[2]);
  k2[2:(T-1)] ~ normal(intercept[2] + slope[2]*k2[1:(T-2)], sigk[2]);
  
  //Model
  snd1 ~ normal(0,1);
  d ~ poisson_log(mpd);
  
}

generated quantities {
  real kf1[F];
  real kf2[F];
  real mmd_f[Co*N*S*F];
  
  kf1[1] = normal_rng( intercept[1] + slope[1]*k1[T-1], sigk[1]);
  kf2[1] = normal_rng( intercept[2] + slope[2]*k2[T-1], sigk[2]);
  for (t in 2:F) kf1[t] = normal_rng(intercept[1] + slope[1]*kf1[t-1], sigk[1]);
  for (t in 2:F) kf2[t] = normal_rng(intercept[2] + slope[2]*kf2[t-1], sigk[2]); 
  for (a in 1:N){
    for (s in 1:S){
      for (t in 1:F){
        for (i in 1:Co){
          mmd_f[Co*(F*(S*(a-1)+s-1)+t-1)+i] = OD_a[i] + G*(s-1) + cOA[i,a]+cDA[i,a]+cO[i]+cD[i] +A_a[a] + A_b[a]*kf1[t] + OD_b[i]*kf2[t];
        }
      }
    }
  }
}
