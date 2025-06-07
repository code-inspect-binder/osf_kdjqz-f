functions {
  real direct_access(int accuracy, real RT, real theta, real P_b, real mu, real delta, real sigma_e){
    
    real log_p_answer_correct;
    real log_p_answer_correct_direct_access;
    real log_p_answer_correct_backtrack;
    real log_p_answer_incorrect;
    
    // theta * (P_b * 1-theta)
    // Log probability of correct answer 
    log_p_answer_correct = log_sum_exp(log(theta), log(P_b) + log1m(theta)); 
    // theta / p_answer_correct
    // Log proportion initial correct answers (without backtracking)
    log_p_answer_correct_direct_access = log(theta) - log_p_answer_correct;
    // (P_b * 1-theta) / p_answer_correct
    // Log Proportion of backtracked answers 
    log_p_answer_correct_backtrack = log(P_b) + log1m(theta) - log_p_answer_correct;
    
    // (1-theta) * (1-P_b)
    // Log probability of incorrect answer
    log_p_answer_incorrect = log1m(P_b) + log1m(theta);
    
    if(accuracy==1) {
      return (log_p_answer_correct + 
            log_sum_exp(
              log_p_answer_correct_direct_access + lognormal_lpdf(RT| mu, sigma_e),
              log_p_answer_correct_backtrack + lognormal_lpdf(RT| mu + delta, sigma_e) ));
    } else {
      return (log_p_answer_incorrect + lognormal_lpdf(RT| mu, sigma_e));
    }
  }

vector direct_access_rng(real theta, real P_b, real mu, real delta, real sigma_e){
    int init_acc;
    int backtrack;
    vector[2] gen;

    init_acc = bernoulli_rng(theta);
    backtrack = 0;
    if (init_acc!=1) backtrack = bernoulli_rng(P_b);
    // Change the answer to 1 if there was backtracking:
    gen[2] = backtrack ? 1 : init_acc;
    { real mu_rng; // it adds the mu_b if there is backtracking:
      mu_rng = mu + (backtrack ? delta : 0);
      gen[1] = lognormal_rng(mu_rng, sigma_e);
    }
    return(gen);
}
}
data {
  int<lower=1> N_obs;                   
  real RT[N_obs];                       
  int<lower=0,upper=1> accuracy[N_obs]; 
  int<lower=-1,upper=1> rctype[N_obs];  
  int<lower=-1,upper=1> group[N_obs];  

  int<lower=1> N_subj;         
  int<lower=1> N_item; 
  int<lower = 1, upper = N_subj> subj[N_obs]; 
  int<lower = 1, upper = N_item> item[N_obs];
  int<lower = 1> n_u; // number of random eff for subj
  int<lower = 1> n_w; // number of random eff for items
   
}

parameters {
  vector[7] beta; //slopes 
  real mu_0;   // intercept 
  real gamma;  //intercept for prob of backtracking in logit space
  real alpha;  // intercept for first initial retrieval in logit space
  
  real<lower=fabs(beta[6])> delta_0;   //effect of backtracking 
  real<lower=fabs(beta[7])> sigma_e_0; //logsd
  
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
  alpha ~ normal(1,.5); //intercept for theta 
  beta ~ normal(0,.5);
  delta_0 ~ normal(0,1); // intercept for the effect of reanalysis
  mu_0 ~ normal(7.5,.6);  // intercept for mu 
  gamma ~ normal(-1,.5);  //intercept for P_b (in logit space)
  sigma_e_0 ~ normal(0,.5);
  tau_u ~ normal(0,.1);
  tau_w ~ normal(0,.1);    
  
  to_vector(z_u) ~  normal(0, .5);
  to_vector(z_w) ~  normal(0, .5);
  L_u ~ lkj_corr_cholesky(2.0);
  L_w ~ lkj_corr_cholesky(2.0);
  
  //log likelihood
  for (i in 1:N_obs){
    //random and fixed effects
    real mu = mu_0 + u[subj[i],1]+ w[item[i],1] + group[i]*beta[1];
    real theta = inv_logit(alpha + u[subj[i],2] + w[item[i],2] + rctype[i]*(beta[2]+ u[subj[i],3]) + group[i]*(beta[3]+ w[item[i],3]) + rctype[i]*group[i]*beta[4]);
    real P_b = inv_logit(gamma + u[subj[i],4] + group[i]*beta[5]); 
    real delta = delta_0 +  group[i]*beta[6];
    real sigma_e = sigma_e_0 + group[i]*beta[7];

    target += direct_access(accuracy[i], RT[i], theta, P_b, mu, delta, sigma_e);
  }
}

generated quantities {
  vector[N_obs] log_lik;
  vector[2] gen;
  vector[N_obs] gen_acc;
  vector[N_obs] gen_rctype;
  vector[N_obs] gen_RT;
  vector[N_obs] gen_group;

  // for first retrieval (theta):
  real mu_i = mu_0 + beta[1];
  real mu_c = mu_0 - beta[1];
  
  real prob_or_i = inv_logit(alpha + beta[2] + beta[3] + beta[4]);
  real prob_or_c = inv_logit(alpha + beta[2] - beta[3] - beta[4]);
  real prob_sr_i = inv_logit(alpha - beta[2] + beta[3] - beta[4]);
  real prob_sr_c = inv_logit(alpha - beta[2] - beta[3] + beta[4]);

  real P_b_i = inv_logit(gamma + beta[5]);
  real P_b_c = inv_logit(gamma - beta[5]);
  
  real delta_i_0 = delta_0 + beta[6];
  real delta_c_0 = delta_0 - beta[6];
  
  real delta_i = exp(mu_i+delta_i_0)-exp(mu_i);
  real delta_c = exp(mu_c+delta_c_0)-exp(mu_c);

  real sigma_e_i = sigma_e_0 + beta[7];
  real sigma_e_c = sigma_e_0 - beta[7];
  

for (i in 1:N_obs){
    //Adjust parameters with random and fixed effects,
    real mu = mu_0 + u[subj[i],1]+ w[item[i],1] + group[i]*beta[1];
    real theta = inv_logit(alpha + u[subj[i],2] + w[item[i],2] + rctype[i]*(beta[2]+ u[subj[i],3]) + group[i]*(beta[3]+ w[item[i],3]) + rctype[i]*group[i]*beta[4]);
    real P_b = inv_logit(gamma + u[subj[i],4] + group[i]*beta[5]); 
    real delta = delta_0 +  group[i]*beta[6];
    real sigma_e = sigma_e_0 + group[i]*beta[7];
    // for loo comparison 
    log_lik[i] = direct_access(accuracy[i], RT[i], theta, P_b, mu, delta, sigma_e);
    // Generate data from the sampled parameters
    gen = direct_access_rng(theta, P_b, mu, delta, sigma_e);
    gen_RT[i] = gen[1];
    gen_acc[i] = gen[2];
    gen_rctype[i] = rctype[i];
    gen_group[i] = group[i];

  }
}
