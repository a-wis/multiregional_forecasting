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
real nu[5];
real sigma[2];
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
}
generated quantities{
  real<lower=0> y[10];
  // real<lower=0> z[2];
  // real<lower=0> sigh;
  
  // sigh = cauchy_rng(0,1);
  for (i in 1:5) {for (j in 1:2) y[(j-1)*5+i]=student_t_rng(nu[i],0,sigma[j]);
  }
  // z[1] = normal_rng(0,sigh);
  // z[2] = student_t_rng(3,0,sigh);
}

