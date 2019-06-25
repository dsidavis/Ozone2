/* 
   Modeling response to Ozone exposure following 
   "Ozone exposure-response model for lung function
   changes: an alternate variability structure" (2013)
   by:
   William F. McDonnell, Paul W. Stewart & Marjo V. Smith

   Model coded in Stan by: Matt Espe
   7 Feb 2019
 */


functions{
  // expand the O3 and Ve to full length
  vector expand_obs(vector x, int[] t){
	int n = max(t);
	vector[n] ans;
	int count = 1;
	for(i in 1:n){
	  ans[i] = x[count];
	  if(i == t[count])
		count += 1;
	  /* print("count: ", count); */
	}
	return ans;
  }

  /* 
	 Find the time pt where the dos threshold is exceeded.
	 After this pt, the formula is significantly simplified
	 This runs on a minute by minute time point, but exits as
	 soon as the change pt has been found. 
  */
  int get_change_pt(vector o3, vector ve, int[] t_stop,
					real e_dos, real thresh){
	int n = max(t_stop);
	vector[n] full_o3 = expand_obs(o3, t_stop);
	vector[n] full_ve = expand_obs(ve, t_stop);
	real cumDr = 0;
	int i = 0;
	
	while(1){
	  i += 1;
	  cumDr += full_o3[i] * full_ve[i] * 1.96;
	  // Once the threshold is exceeded, i = the start pt
	  if(1/(1 + exp(-20 * (i - (e_dos * (i/cumDr))))) > thresh)
		break;
	  
	  // if the threshold is never exceeded, full zeros?
	  if(i == n)
		break;
	  
	}
	
	return i;
  }

  /* 
	 Function to get change in intermediate variable X by time
   */
  vector get_X(vector O3, vector Ve,
				 int[] Time,
				 real e_dos, real e_k){
	int N = dims(O3)[1];
	vector[N] X = rep_vector(0.0, N);
	real Vm;
	real Ta;
	real Tb;
	real TD;
	real X_previous;
	int start_pt = get_change_pt(O3, Ve, Time, e_dos, 0.001); 
	/* print("change point is: ", start_pt); */
	
	for(j in 1:N){
	  // if the cumDR is less than DOS threshold, X must be 0
	  if(Time[j] < start_pt){
		if(j < N){
		  if(Time[j+1] > start_pt
		}
		/* print("Under start pt"); */
		X[j] = 0;
	  } else {
		real DR = O3[j] * Ve[j] * 1.96;
		// Once the DOS threshold is passed, use formula
		Ta = j > 1 ? Time[j-1] : 0.0;
		Tb = Time[j];
		TD = (Tb - Ta);
		X_previous = j > 1 ? X[j-1] : 0;
		
		X[j] = X_previous *
		  exp(e_k * TD) +
		  (DR * (1 - exp(-e_k * TD))) / e_k;
	  }
	}
	/* print("X = ", X); */
	return X;
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
  // These are padded with zeros
  vector[max_timepts] Ve[n_obs];
  vector[max_timepts] Cm[n_obs];
  int Time[n_obs, max_timepts];
  int dFEV1_measure_idx[n_obs, max_n_dFEV1];
  vector[max_n_dFEV1] dFEV1[n_obs];
  real dos;
  real k;
  real a;
  real sigma;
}

transformed data{
  real e_dos = exp(dos);
  real e_k = exp(k);
  real e_a = exp(a);
}

parameters{
/*   real<lower = 1.5, upper = 8> log_dos; */
/*   real<lower = -20, upper = -1> log_k; */
/*   real<lower= -20, upper = 2> log_a; */
/*   real<lower = 0> sigma; */
}

model{
/*   for(n in 1:n_obs){ */
/* 	int idx = n_timepts[n]; */
/* 	vector[idx] X = get_X(Cm[n][:idx],  */
/* 						  Ve[n][:idx], */
/* 						  Time[n][:idx], */
/* 						  dos, k); */
	
/* 	int comp_idx[n_dFEV1[n]] = dFEV1_measure_idx[n][:n_dFEV1[n]]; */

/* 	/\* print("Predicted: ", X[comp_idx] * exp(log_a)); *\/ */
/* 	/\* print("Observed: ", dFEV1[n][:n_dFEV1[n]]); *\/ */
/* 	// Likelihood */
/* 	dFEV1[n][:n_dFEV1[n]] ~ normal(X[comp_idx] * exp(a), sigma); */
/*   } */
}

generated quantities{
  real log_lik = 0;
  real aic;

  for(n in 1:n_obs){
	int idx = n_timepts[n];
	vector[idx] X = get_X(Cm[n][:idx], 
						  Ve[n][:idx],
						  Time[n][:idx],
						  e_dos, e_k);
	
	int comp_idx[n_dFEV1[n]] = dFEV1_measure_idx[n][:n_dFEV1[n]];

  	// Likelihood
  	log_lik += normal_lpdf(dFEV1[n][:n_dFEV1[n]] | X[comp_idx] * e_a, sigma);
  }
  
  // k = number of betas + number of random effects
  aic = 2 * 3 - 2 * log_lik;
}
