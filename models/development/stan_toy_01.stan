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
  matrix [R-1, R-1] mat_Rc;
  matrix [N-1, N-1] mat_Ra;
  real sig_t;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[N-1] ben1;
  vector[R-1] ben2;
  real<lower=0> sig1;
  real<lower=0> sig2[2];
  real<lower=0> sigk[2];
  vector[T-1] k1; //time series
  vector[T-1] k2; //time series 
  // real<lower=0, upper=1> slope[2]; 
  // real intercept[2];
  // real snd1[R*N*T];
  real A_arnd[N];
  real R_arnd[R];
}

transformed parameters {
  vector[N] A_a;
  vector[N] A_b;
  vector[R] R_b;
  vector[R] R_a;
  real mmd[R*N*T]; 

  
  for (i in 1:N) A_a[i] = A_arnd[i] * sig2[1];
  for (r in 1:R) R_a[r] = R_arnd[r] * sig2[2];
  
  for (i in 1:(N-1)){  
    A_b[i] = ben1[i];
  }
  A_b[N] = 1-sum(ben1[1:(N-1)]);
  
  for (i in 1:(R-1)){  
    R_b[i] = ben2[i];
  }
  R_b[R] = 1-sum(ben2[1:(R-1)]);
  
  for (a in 1:N){
    for (t in 1:T){
      for (i in 1:R){
        mmd[R*(T*(a-1)+t-1)+i] = R_a[i] + A_a[a] + (t == 1 ? 0 : A_b[a]*k1[t-1]) + (t == 1 ? 0 : R_b[i]*k2[t-1]);
        
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
  vector[N-1] mb1;
  vector[R-1] mb2;
  matrix [N-1, N-1] mat_R1;
  matrix [R-1, R-1] mat_R2;
  real taub1;
  real taub2;
  
  // priors for precision
  taub1 = N*N;
  taub2 = R*R;
  sig1  ~ student_t(2.5,0,sig_t);
  sigk[1]  ~ student_t(2.5,0,sig_t);
  sigk[2]  ~ student_t(2.5,0,sig_t);
  sig2  ~ student_t(2.5,0,sig_t);

  // priors for random effects  

  // prior for A_b  
  for (i in 1:(N-1)) mb1[i] = 1.0/N;
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      mat_R1[i,j] = mat_Ra[i,j]*taub1;
    }
  }
  ben1 ~ multi_normal_prec(mb1,mat_R1);
  
  // prior for OD_b  
  for (i in 1:(R-1)) mb2[i] = 1.0/R;
  for (i in 1:(R-1)){
    for (j in 1:(R-1)){
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
  // snd1 ~ normal(0,1);
  // Grnd ~ normal(0,1);
  A_arnd ~ normal(0,1);
  R_arnd ~ normal(0,1);
  
  y ~ normal(mmd, sig1);
}

