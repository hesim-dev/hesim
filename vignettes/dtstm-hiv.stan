functions {
 // Function to apply a relative risk to a transition probability matrix
 matrix adjust_p(matrix p, real rr){
    int n_cols = cols(p);
    int n_rows = rows(p);
    matrix[n_rows, n_cols] p_new = p;
    
    // Adjust each row of matrix using relative risks
    for (i in 1:n_rows){
       p_new[i, (i + 1):n_cols] = rr * p[i, (i + 1):n_cols];
       p_new[i, i] = 1 - sum(p_new[i, (i + 1):n_cols]);
    }
    return(p_new);
 }
}

data {
  int<lower=0> n_strategies; // Number of treatment strategies
  int<lower=0> n_patients; // Number of patients
  int<lower=0> n_times; // Number of time intervals
  matrix[4, 4] alpha; // Concentration parameters of Dirichlet distributions
  real lrr_mean; // Mean of the log relative risk
  real lrr_se; // Standard error of the log relative risk
}

transformed data{
  matrix[4, 4] t_alpha = alpha'; // Transpose alpha. This is needed to randomly
                                // sample from the Dirichlet distribution in Stan.
}

 generated quantities {
   matrix[4, 4] p_mono; // Transition probability matrix for monotherapy
   matrix[4, 4] p_combo; // Transition probability matrix for combination therapy
   matrix[4, 4] p[n_strategies, n_patients, n_times]; // 4D array to pass to hesim
   
   // Draw relative risk
   real rr = lognormal_rng(lrr_mean, lrr_se); 
   
   // Draw baseline transition probability matrix for monotherapy
   for (i in 1:4){
     p_mono[i, ] = dirichlet_rng(col(t_alpha, i))'; 
   }                                                
   
   // Compute transition probability matrix for combination therapy using 
   // relative risk
   p_combo = adjust_p(p_mono, rr);

   // Store transition probabilitiy matrices in 4D array 
   p[1, 1, 1] = p_mono;  // Monotherapy, time interval 1 
   p[2, 1, 1] = p_combo; // Combination therapy, time interval 1
   p[1, 1, 2] = p_mono;  // Monotherapy, time interval 2
   p[2, 1, 2] = p_mono;  // Combination therapy, time interval 2 (treatment is
                         // assumed to only last 2 years) 
}
