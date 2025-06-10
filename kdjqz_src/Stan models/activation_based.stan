functions {
  real race(int winner, real RT, real accum_SR_int_mu, real accum_OR_int_mu, real accum_sigma){
    
    real log_lik;
    log_lik = 0;
    
    if(winner==1){ // SR are coded as 1 in "winner"
      log_lik += lognormal_lpdf(RT| accum_SR_int_mu, accum_sigma); 
      log_lik += lognormal_lccdf(RT| accum_OR_int_mu, accum_sigma); 
    }
    else {
      log_lik += lognormal_lpdf(RT| accum_OR_int_mu, accum_sigma);
      log_lik += lognormal_lccdf(RT| accum_SR_int_mu, accum_sigma); 
    }
    return(log_lik);
  }
  
  // RTs for ppc
  vector race_rng(real mu_1, real mu_2, real sigma, int rctype){
    vector[2] gen;
    real accum_SR_RT = lognormal_rng(mu_1, sigma);
    real accum_OR_RT = lognormal_rng(mu_2, sigma);
    
    if(accum_SR_RT < accum_OR_RT){
      gen[1] = accum_SR_RT;
      if(rctype == -1){ //SRs are coded as -1
        gen[2] = 1;
      }
      else {
        gen[2] = 0;
      }
    }
    else {
      gen[1] = accum_OR_RT;
      if(rctype == 1){ //ORs are coded as 1
        gen[2] = 1;
      }
      else {
        gen[2] = 0;
      }
    }
    return(gen);
  }
}
data { 
  int<lower = 1> N_obs; 
  int<lower = 1> N_choices; 
  int<lower = 1> n_u; //number of random effects for subj
  int<lower = 1> n_w; // number of random effects for items
  int<lower =-1, upper = 1> rctype[N_obs];
  int<lower =-1, upper = 1> group[N_obs];
  int<lower = 1, upper = N_choices> winner[N_obs];
  vector<lower = 0>[N_obs] RT;
  
  int<lower = 1> N_subj;
  int<lower = 1> N_item;
  int<lower=1> subj[N_obs];    
  int<lower=1> item[N_obs]; 
}
parameters{
  vector[7] beta;
  real alpha[N_choices];
  real<lower=fabs(beta[7])> sigma_e;
  
  cholesky_factor_corr[n_u] L_u;  
  cholesky_factor_corr[n_w] L_w;  
  vector<lower=0>[n_u] tau_u; 
  vector<lower=0>[n_w] tau_w;
  matrix[n_u,N_subj] z_u;           
  matrix[n_w,N_item] z_w;     
}
transformed parameters {
  matrix[N_subj, n_u] u;
  matrix[N_item, n_w] w;
  u = (diag_pre_multiply(tau_u, L_u) * z_u)';
  w = (diag_pre_multiply(tau_w, L_w) * z_w)';
}
  
model {
  //priors
  target += normal_lpdf(alpha| 7.5,.6);
  target += normal_lpdf(beta |    0 ,.5);
  target += normal_lpdf(sigma_e | 0 ,.5)- normal_lccdf(0 | 0,.5);
  target += normal_lpdf(tau_u | 0 ,.1)  - normal_lccdf(0 | 0,.1);
  target += normal_lpdf(tau_w | 0 ,.1)  - normal_lccdf(0 | 0,.1);
  
  target += lkj_corr_cholesky_lpdf(L_u | 2);
  target += lkj_corr_cholesky_lpdf(L_w | 2);
  target += normal_lpdf(to_vector(z_u) | 0, .5);
  target += normal_lpdf(to_vector(z_w) | 0, .5);

  for (n in 1:N_obs) {
    real accum_SR_int_mu = alpha[1] + u[subj[n],1] + w[item[n],1] + group[n]*(beta[1]+w[item[n],3]) + rctype[n]*(beta[3]+u[subj[n],3]) + group[n]*rctype[n]*beta[5];
    real accum_OR_int_mu = alpha[2] + u[subj[n],2] + w[item[n],2] + group[n]*(beta[2]+w[item[n],4]) + rctype[n]*(beta[4]+u[subj[n],4]) + group[n]*rctype[n]*beta[6];
    real accum_sigma = sigma_e + group[n]*beta[7];
        target += race(winner[n], RT[n], accum_SR_int_mu, accum_OR_int_mu, accum_sigma);
  }
}

generated quantities {
  vector[N_obs] rt_1;
	vector[N_obs] rt_2;
  vector[N_obs] gen_acc;
  vector[N_obs] gen_rctype;
  vector[N_obs] gen_RT;
  vector[N_obs] gen_subj;
  vector[N_obs] gen_group;
  vector[N_obs] log_lik;
  
  real sigma_i = sigma_e + beta[7];
  real sigma_c = sigma_e - beta[7];
  
  for (n in 1:N_obs){
    vector[2] gen;  
    real mu_1;
    real mu_2;
    real sigma;
    
    mu_1 = alpha[1] + u[subj[n],1] + w[item[n],1] + group[n]*(beta[1]+w[item[n],3]) + rctype[n]*(beta[3]+u[subj[n],3]) + group[n]*rctype[n]*beta[5];
    mu_2 = alpha[2] + u[subj[n],2] + w[item[n],2] + group[n]*(beta[2]+w[item[n],4]) + rctype[n]*(beta[4]+u[subj[n],4]) + group[n]*rctype[n]*beta[6];
    sigma = sigma_e + group[n]*beta[7]; 
    
    gen = race_rng(mu_1, mu_2, sigma, rctype[n]);
    gen_RT[n] = gen[1]; // estimated listening times 
    gen_acc[n] = gen[2];
    gen_rctype[n] = rctype[n];
    gen_subj[n] = subj[n];
    gen_group[n] = group[n];
    
    rt_1[n] = lognormal_rng(mu_1, sigma); // finishing times accumulator SR
		rt_2[n] = lognormal_rng(mu_2, sigma); // finishing times accumulator OR
		
    log_lik[n] = race(winner[n], RT[n], mu_1, mu_2, sigma);
  }
} 
