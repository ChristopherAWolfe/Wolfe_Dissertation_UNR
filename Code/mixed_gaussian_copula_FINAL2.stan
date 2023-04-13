functions{
  
  // All associated code is adapted from: 
  // https://spinkney.github.io/helpful_stan_functions/group__gaussian__copula.html
  
  // There are 2 functions that define a Gaussian Copula in Stan. First, `multi_
  // normal_cholesky_copula_lpdf` and second, `centered_gaussian_copula_choelsky_`
  // These work in tandom to first aggregate all the marginal calculations and 
  // jacobian adjustment and then increment the log-probability.
  
  real multi_normal_cholesky_copula_lpdf(matrix U, matrix L) {
    int N = rows(U);
    int J = cols(U);
    matrix[J, J] Gammainv = chol2inv(L);
    return -N * sum(log(diagonal(L))) - 0.5 * sum(add_diag(Gammainv, -1.0) .* crossprod(U));
  }
  
  real centered_gaussian_copula_cholesky_(array[,] matrix marginals, matrix L) {
    // Extract dimensions of final outcome matrix
    int N = rows(marginals[1][1]);
    int J = rows(L);
    matrix[N, J] U;
  
    // Iterate through marginal arrays, concatenating the outcome matrices by column
    // and aggregating the log-likelihoods (from continuous marginals) and jacobian
    // adjustments (from discrete marginals)
    real adj = 0;
    int pos = 1;
    for (m in 1 : size(marginals)) {
      int curr_cols = cols(marginals[m][1]);
    
      U[ : , pos : (pos + curr_cols - 1)] = marginals[m][1];
    
      adj += sum(marginals[m][2]);
      pos += curr_cols;
    }
  
    // Return the sum of the log-probability for copula outcomes and jacobian adjustments
    return multi_normal_cholesky_copula_lpdf(U | L) + adj;
  }

  // There are 2 function types that define each univariate marginal distribution depending 
  // on data type: polychotomous ordinal or continuous.The ordinal is adapted into
  // two separate functions depending on data structure: (N x K matrix vs. N vector)
  // The ordinal data follows the data augmentation approach highlighted in Albert
  // & Chib, 1993. The continuous marginal is std_normal.
  
  array[] matrix ordered_probit_marginal(array[,] int y, matrix mu_glm, matrix u_raw, array[] vector cutpoints) {
    int N = rows(mu_glm); // number of observations
    int K = cols(mu_glm); // number of polytomous variables
    array[2] matrix[N, K] rtn; // empty 2D array to return

    for(k in 1:K){
      for(n in 1:N){
        int C = num_elements(cutpoints[k,]) + 1; // total number of ord categories
        if(y[n,k] == 99){ // missing data
          rtn[1][n,k] = inv_Phi(u_raw[n,k]); // missing RV
          rtn[2][n,k] = 0;
        } else if(y[n,k] == 1){ // lowest bound
          real bound = Phi((cutpoints[k,1] - mu_glm[n,k])); // data augmentation
          rtn[1][n,k] = inv_Phi((bound*u_raw[n,k])); // latent RV
          rtn[2][n,k] = log(bound); // jacobian
        } else if (y[n,k] == C){ // highest bound
          real bound = Phi((cutpoints[k, C - 1] - mu_glm[n,k])); // data augmentation
          rtn[1][n,k] = inv_Phi(bound + (1-bound)*u_raw[n,k]); // latent RV
          rtn[2][n,k] = log1m(bound); // jacobian
        } else { // in between 
          real ub = Phi((cutpoints[k ,y[n,k]] - mu_glm[n,k])); // data augmentation
          real lb = Phi((cutpoints[k, y[n,k] - 1] - mu_glm[n,k])); // data augmentation
          rtn[1][n,k] = inv_Phi((lb + (ub-lb)*u_raw[n,k])); // latent RV
          rtn[2][n,k] = log(ub-lb); // jacobian
        }
      }
    }
  return rtn;
  }
  
  array[] matrix ordered_probit_marginal_uni(int[] y, real[] mu_glm, real[] u_raw, vector cutpoints) {
    int N = num_elements(mu_glm); // number of observations
    array[2] matrix[N, 1] rtn; // empty 2D array to return
    
    for(n in 1:N){
      int C = num_elements(cutpoints) + 1; // total number of ord categories
      if(y[n] == 99){ // missing data
        rtn[1][n,1] = inv_Phi(u_raw[n]); // missing RV
        rtn[2][n,1] = 0;
        } else if(y[n] == 1){ // lowest bound
        real bound = Phi((cutpoints[1] - mu_glm[n])); // data augmentation
        rtn[1][n,1] = inv_Phi((bound*u_raw[n])); // latent RV
        rtn[2][n,1] = log(bound); // jacobian
      } else if (y[n] == C){ // highest bound
        real bound = Phi((cutpoints[C - 1] - mu_glm[n])); // data augmentation
        rtn[1][n,1] = inv_Phi(bound + (1-bound)*u_raw[n]); // latent RV
        rtn[2][n,1] = log1m(bound); // jacobian
      } else { // in between 
        real ub = Phi((cutpoints[y[n]] - mu_glm[n])); // data augmentation
        real lb = Phi((cutpoints[y[n] - 1] - mu_glm[n])); // data augmentation
        rtn[1][n,1] = inv_Phi((lb + (ub-lb)*u_raw[n])); // latent RV
        rtn[2][n,1] = log(ub-lb); // jacobian
      }
    }
    return rtn;
  }

  matrix[] normal_marginal(matrix y, matrix mu_glm, matrix sigma) {
    int N = rows(mu_glm); // number of observations
    int J = cols(mu_glm); // number of continuous variables
    matrix[N, J] rtn[2]; // empty 2D array to return
    // Initialise the jacobian adjustments to zero, as vectorised lpdf will be used
    rtn[2] = rep_matrix(0, N, J);

    for (j in 1 : J) {
      rtn[1][ : , j] = (y[ : , j] - mu_glm[ : , j]) ./ sigma[ :,j]; // center RV
      rtn[2][1, j] = normal_lpdf(y[ : , j] | mu_glm[ : , j], sigma[ :,j]); // "jacobian"
    }
  return rtn;
  }

}
data{
  
  // size variables //
  int N; // total # of observations (rows)
  int J; // total # of continuous observations
  int K_dentition; // total # of dental variables
  int K_ef; // total # of ef variables 
  int K_pelvis; // total # of pelvic variables
  int K_all; // all ordinal variables
  
  //continuous prep //
  int present; // total # of complete data in continuous responses
  int missing; // total # of missing data in continuous responses
  real y_continuous[present]; // vector of reals of responses across all vars
  int index_pres[present, 2]; // Matrix (row, col) of non-missing value indices
  int index_miss[missing, 2]; // Matrix (row, col) of missing value indices
  
  // polychotomous prep //
  int y_dentition[N,K_dentition]; // 2D array of integers n obs by k_dentition
  int y_ef[N,K_ef]; // 2D array of integers n obs by k_ef
  int y_pelvis[N,K_pelvis]; // 2D array of integers n obs by k_pelvis
  int y_carpal[N]; // 2D array of integers n obs by k_carpal
  int y_tarsal[N]; // 2D array of integers n obs by k_tarsal
  
  // predictor (age) //
  real X[N]; //scaled

}
parameters{
  
  // continuous parameters //
  real y_missing_cont[missing]; // missing data parameter the size of all missing vars
  vector<lower=0>[J] constant; //constant in continuous mean function
  vector<lower=0>[J] exponent; // exponent in continuous mean function
  vector[J] offset; // offset in continuous mean function
  vector<lower=0>[J] noise1; // param of linear noise function
  vector<lower=0>[J] noise2; // param of linear noise function

  // Polychotomous parameters //
  vector<lower=0>[K_dentition] constant_dentition;
  vector<lower=0>[K_ef] constant_ef;
  vector<lower=0>[K_pelvis] constant_pelvis;
  real<lower=0> constant_carpal;
  real<lower=0> constant_tarsal;
  matrix<lower=0, upper=1>[N,K_dentition] u_dentition; // nuisance dentition
  matrix<lower=0, upper=1>[N,K_ef] u_ef; // nuisance ef
  matrix<lower=0, upper=1>[N,K_pelvis] u_pelvis; // nuisance pelvis
  real<lower=0, upper=1> u_carpal[N]; // nuisance carpal
  real<lower=0, upper=1> u_tarsal[N]; // nuisance tarsal
  ordered[11] threshold_dentition[16]; // dental thresholds
  ordered[6] threshold_ef[16]; // ef thresholds
  ordered[2] threshold_pelvis[2]; // pelvis thresholds
  ordered[8] threshold_carpal; // carpal thresholds
  ordered[7] threshold_tarsal; // tarsal thresholds
  
  // Cholesky factor (copula parameter) //
  cholesky_factor_corr[J+K_all] L;
  
}
transformed parameters{
  
  // Parameter Declarations //
  matrix[N,J] mu_continuous; // continuous mean
  matrix[N,J] sd_continuous; // continuous sd
  matrix[N,J] y_full_continuous; // full y matrix including missing data parameters
  matrix[N,K_dentition] mu_dentition; // dental mean
  matrix[N,K_ef] mu_ef; // ef mean
  matrix[N,K_pelvis] mu_pelvis; // pelvis mean
  real mu_carpal[N]; // carpal mean
  real mu_tarsal[N]; // tarsal mean

  // Continuous Parameters //
  for(i in 1:N){
    for(j in 1:J){
      mu_continuous[i,j] = constant[j]*X[i]^exponent[j] + offset[j];
      sd_continuous[i,j] = noise1[j]*(1 + X[i]*noise2[j]);
    }
  }
  // Fill y_full continuous with present data
  for(n in 1:present) {
      y_full_continuous[index_pres[n,1]][index_pres[n,2]] = y_continuous[n];
    }
  // Fill the rest of y_full continuous with missing value "parameters"
  for(n in 1:missing){
      y_full_continuous[index_miss[n,1]][index_miss[n,2]] = y_missing_cont[n];
    }

 // Polytomous Parameters //
 for(i in 1:N){
   for(k in 1:K_dentition){
     mu_dentition[i,k] = constant_dentition[k]*X[i];
    }
  }
 for(i in 1:N){
   for(k in 1:K_ef){
     mu_ef[i,k] = constant_ef[k]*X[i];
    }
  }
 for(i in 1:N){
   for(k in 1:K_pelvis){
     mu_pelvis[i,k] = constant_pelvis[k]*X[i];
    }
  }
 for(i in 1:N){
     mu_carpal[i] = constant_carpal*X[i];
  }
 for(i in 1:N){
     mu_tarsal[i] = constant_tarsal*X[i];
  }
  
}
model{
  
  // Priors // 
  constant ~ std_normal();
  exponent ~ std_normal();
  offset ~ std_normal();
  noise1 ~ std_normal();
  noise2 ~ std_normal();
  constant_dentition ~ std_normal();
  constant_ef ~ std_normal();
  constant_pelvis ~ std_normal();
  constant_carpal ~ std_normal();
  constant_tarsal ~ std_normal();

  // Threshold Priors //
  for(i in 1:K_dentition){
    for(t in 1:11){
      threshold_dentition[i,t] ~ normal(t+1,1);
    }
  }
  for(i in 1:K_ef){
    for(t in 1:6){
      threshold_ef[i,t] ~ normal(t+1,1);
    }
  }
  for(i in 1:K_pelvis){
    for(t in 1:2){
      threshold_pelvis[i,t] ~ normal(t+1,1);
    }
  }
  
    for(t in 1:8){
      threshold_carpal[t] ~ normal(t+1,1);
    }
  
    for(t in 1:7){
      threshold_tarsal[t] ~ normal(t+1,1);
    }

  
  
  // Copula Parameter Prior //
  L ~ lkj_corr_cholesky(8);
  
  // Likelihood //
  target += centered_gaussian_copula_cholesky_(
    {normal_marginal(y_full_continuous, mu_continuous, sd_continuous),
    ordered_probit_marginal(y_dentition, mu_dentition, u_dentition, threshold_dentition),
    ordered_probit_marginal(y_ef, mu_ef, u_ef, threshold_ef),
    ordered_probit_marginal(y_pelvis, mu_pelvis, u_pelvis, threshold_pelvis),
    ordered_probit_marginal_uni(y_carpal, mu_carpal, u_carpal, threshold_carpal),
    ordered_probit_marginal_uni(y_tarsal, mu_tarsal, u_tarsal, threshold_tarsal)}, L);
  
}
generated quantities{
  
  // Here I put the correlation matrix back together for ease of interpretation
  corr_matrix[J+K_all] cor_mat = multiply_lower_tri_self_transpose(L);
}
