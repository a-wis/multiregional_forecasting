//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> R;         // origins and destinations (R x R-1)
  int<lower=0> N;         // Age
  int<lower=0> T;         // Time
  vector[R*N*T] y;
  vector[N] A_b;
  vector[R] R_b;
  real sig_t;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[N] A_a;
  vector[R] R_a;
  real<lower=0> sig1;
  // real<lower=0> sig2[2];
  real<lower=0> sigk[2];
  real<lower=0, upper=1> slope[2];
  // vector[T] kn1; //time series
  // vector[T] kn2; //time series 
  real intercept[2];
  // real snd1[R*N*T];
  // real A_arnd[N];
  // real R_arnd[R];
  real k1_rnd[T];
  real k2_rnd[T];
}

transformed parameters {
  // vector[N] A_a;
  // vector[R] R_a;
  vector[T] k1; //time series
  vector[T] k2; //time series 
  real mmd[R*N*T]; 
  // real mpd[R*N*T];
  
  // for (i in 1:(T-1)) k1[i] = kn1[i] - ;

  // for (i in 1:N) A_a[i] = A_arnd[i] * sig2[1];
  // for (r in 1:R) R_a[r] = R_arnd[r] * sig2[2];

  k1[1] = intercept[1]+k1_rnd[1]*sigk[1];
  k2[1] = intercept[2]+k2_rnd[1]*sigk[2];
  for (i in 2:(T)) {
    k1[i] = intercept[1]+slope[1]*k1[i-1] + k1_rnd[i]*sigk[1];
    k2[i] = intercept[2]+slope[2]*k2[i-1] + k2_rnd[i]*sigk[2];
  }
  
  
  for (a in 1:N){
    for (t in 1:T){
      for (i in 1:R){
        mmd[R*(T*(a-1)+t-1)+i] = R_a[i] + A_a[a] + A_b[a]*k1[t] + R_b[i]*k2[t];
        
      }
    }
  }

  
  // for (i in 1:(R*N*T)){
  //   mpd[i] = mmd[i] + sig1*snd1[i];
  // }  
  
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // priors for precision
  sig1  ~ student_t(2.5,0,sig_t);
  sigk[1]  ~ student_t(2.5,0,sig_t);
  sigk[2]  ~ student_t(2.5,0,sig_t);
  // sig2  ~ student_t(2.5,0,sig_t);

  // priors for random effects  

  // prior for A_b  

  A_a ~ normal(0,5);
  R_a ~ normal(0,5);

  // time series model AR1
  // intercept ~ normal(0,10);
  // slope ~ normal(0,1);
  // 
  // k1[1] ~ normal( 0.0, sigk[1]);
  // k1[2:(T-1)] ~ normal(k1[1:(T-2)], sigk[1]);
  // k2[1] ~ normal( 0.0, sigk[2]);
  // k2[2:(T-1)] ~ normal(k2[1:(T-2)], sigk[2]);
  
  //Model
  // snd1 ~ normal(0,1);
  // Grnd ~ normal(0,1);
  // A_arnd ~ normal(0,1);
  // R_arnd ~ normal(0,1);
  k1_rnd ~ normal(0,1);
  k2_rnd ~ normal(0,1);
  y ~ normal(mmd, sig1);
}
// generated quantities{
//   real km1;
//   real km2;
//   
//   km1=mean(k1[1:T]);
//   km2=mean(k2[1:T]);
// }
