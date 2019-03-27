/* Modeling response to Ozone exposure following 
   "Ozone exposure-response model for lung function
   changes: an alternate variability structure" (2013)
   by:
   William F. McDonnell, Paul W. Stewart & Marjo V. Smith

   Model coded in Stan by: Matt Espe
   7 Feb 2019
 */


functions{
  /* 
	 Function to get change in intermediate variable X by time
   */
  vector get_XB5(vector Cm, vector Cs, real Vs, vector Ve,
				 real BSA, vector Time,
				 real B5, real B6){
	int N = dims(Cm)[1];
	vector[N] XB5 = rep_vector(0.0, N);
	real Vm;
	real Ta;
	real Tb;
	real TD;
	real XB5_previous;
	
	for(j in 1:N){
	  Vm = (Ve[j] ^ B6) / (BSA ^ B6);
	  Ta = j > 1 ? Time[j-1] : 0.0;
	  Tb = Time[j];
	  TD = (Tb - Ta);
	  XB5_previous = j > 1 ? XB5[j-1] : 0;
	  
	  XB5[j] = XB5_previous * (exp(-B5 * TD)) +
		(Cm[j] * Vm * (B5^-1)) * (1 - exp(-B5 * TD)) +
		(Cm[j] * Vs * (B5^-2)) * (((1 - B5 * Ta) * exp(-B5 * TD)) - (1 -B5 * Tb)) +
		(Cs[j] * Vm * (B5^-2)) * (((1 - B5 * Ta) * exp(-B5 * TD)) - (1 -B5 * Tb)) +
		(Cs[j] * Vs * (B5^-3)) * (((-2 + (2 * B5 * Ta) -
								 (B5^2 * Ta^2)) * exp(-B5 * TD)) -
							   (-2 + (2 * B5 * Tb) - (B5^2 * Tb^2)));
    }
	return XB5;
  }

  /*
	Takes change in intermediate variable X, and adds adjustments for 
	age and BMI, each centered to their mean values.
   */
  vector get_pop_median (vector XB5, real age_c, real BMI_c, 
						 real B1, real B2, real B3,
						 real B4, real B8, real B9) {
	int N = dims(XB5)[1];
	vector[N] XB5G;
	real F1 = B1 + B2 * age_c + B8 * BMI_c;
	vector[N] T1;
	real T2 = 1 + B4;
	vector[N] Median;
	
	for(n in 1:N)
	  XB5G[n] = XB5[n] >= B9 ? XB5[n] - B9 : 0;
	T1 = 1 + B4 * exp(-B3 * XB5G);
	for(n in 1:N)
	  Median[n] = F1 * (1/T1[n] - 1/T2);
	
	return Median;
  }
}

data{
  int max_timepts; // max number of time points
  int max_n_dFEV1;
  int n_obs; // individual x study x exposure = obs
  int n_ind; // individuals
  int n_dFEV1[n_obs]; // n dFEV1 measurements per obs
  int n_timepts[n_obs]; // n timepts per obs
  int ind[n_obs]; 
  vector[n_obs] age;
  vector[n_obs] BMI;
  vector[n_obs] BSA;
  // These are padded with zeros
  vector[max_timepts] Ve[n_obs];
  vector[max_timepts] Cm[n_obs];
  vector[max_timepts] Cs[n_obs];
  vector[max_timepts] Time[n_obs];
  int dFEV1_measure_idx[n_obs, max_n_dFEV1];
  vector[max_n_dFEV1] dFEV1[n_obs];
  real<lower = 0> sigma_U; // prior sd of random effects
}

transformed data{
  // Center each according to values in paper
  real Vs = 0; // Confused about why this is hardcoded to == 0
  vector[n_obs] age_c = age - 23.8;
  vector[n_obs] BMI_c = BMI - 23.1;

}

parameters{
  vector[9] B;
  vector[n_ind] U; // random effects by ind.
  real sigma;
}

model{
  // distribution of ind. random effects
  U ~ normal(0, sigma_U);
  B ~ normal(0, 100);
  
  for(n in 1:n_obs){
	int idx = n_timepts[n];
	vector[idx] XB5 = get_XB5(Cm[n][:idx], Cs[n][:idx],
									   Vs, Ve[n][:idx],
									   BSA[n], Time[n][:idx],
									   B[5], B[6]);
	vector[idx] med = get_pop_median(XB5, age_c[n], BMI_c[n],
									 B[1], B[2], B[3], B[4], B[8], B[9]);

	int comp_idx[n_dFEV1[n]] = dFEV1_measure_idx[n][:n_dFEV1[n]];

	// Likelihood
	dFEV1[n][:n_dFEV1[n]] ~ normal(med[comp_idx] * exp(U[ind[n]]), sigma);
  }
}

generated quantities{
  real log_lik = 0;
  real aic;

  for(n in 1:n_obs){
	int idx = n_timepts[n];
	vector[idx] XB5 = get_XB5(Cm[n][:idx], Cs[n][:idx],
									   Vs, Ve[n][:idx],
									   BSA[n], Time[n][:idx],
									   B[5], B[6]);
	vector[idx] med = get_pop_median(XB5, age_c[n], BMI_c[n],
									 B[1], B[2], B[3], B[4], B[8], B[9]);

	int comp_idx[n_dFEV1[n]] = dFEV1_measure_idx[n][:n_dFEV1[n]];

	// Likelihood
	log_lik += normal_lpdf(dFEV1[n][:n_dFEV1[n]] | med[comp_idx] * exp(U[ind[n]]), sigma);
  }
  
  // k = number of betas + number of random effects
  aic = 2 * (8 + n_ind) - 2 * log_lik; 
}
