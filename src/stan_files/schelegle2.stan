/*
  Model for change in %FEV1 after exposure to O3

  Matt Espe
  Aug 2018
*/

functions{
  
  vector seq(int start, int end){
	int N = (1 + end) - start;
	vector[N] x;
	for(i in 1:N)
	  x[i] = start + i - 1;
	return x;	
  }
  
  real UOS(real DR, real CumFrDos, real Dos, real t) {
	real UOS = DR / (1 + exp(-20 * (t - (t / CumFrDos))));
	
	return UOS;
  }
  
  real deltaX(real uos, real k, real x_prev){
	real r = (1 - exp(-k));
	real ceoss = uos / k;
	real x = x_prev + ((ceoss - x_prev) * r);
	
	return x;
  }

  int int_max(vector x){
	int ans = 0;
	while(ans < max(x))
	  ans += 1;
	return ans;
  }
  
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
  
  vector experimentFEV1(vector o3, vector ve, int[] t_stop,
						real dos, real k, real a){
	int n = max(t_stop);
	vector[n] full_o3 = expand_obs(o3, t_stop);
	vector[n] full_ve = expand_obs(ve, t_stop);
	
	vector[n] dFEV1 = rep_vector(0, n);
	vector[n] x = rep_vector(0, n);
	real dr = 0;
	real CumFrDos = 0;
	real FrDos = 0;

	/* print("full o3: ", full_o3); */
	for(i in 2:n){
	  real uos;
	  dr = full_o3[i] * full_ve[i] * 1.96;
	  FrDos = dr / dos;
	  CumFrDos += FrDos;
	  /* print("CumFrDos: ", CumFrDos); */
	  uos = UOS(dr, CumFrDos, dos, round(i));
	  /* print("UOS: ", uos); */
	  x[i] = deltaX(uos, k, x[i-1]);
	  print("x = ", x[i]);
	  print("t = ", i);
	}
	dFEV1 = x * a;
	/* print("dFEV1 = ", dFEV1); */
	return dFEV1;
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
  vector[max_timepts] Cs[n_obs];
  int Time[n_obs,max_timepts];
  int dFEV1_measure_idx[n_obs, max_n_dFEV1];
  vector[max_n_dFEV1] dFEV1[n_obs];
  real k;
  real dos;
  real a;
  /* real<lower = 0> sigma; */
}
transformed data{
}
parameters{
  /* real<lower = 1.6, upper = 7.9> dos; */
  /* real<lower = -23, upper = -1.7> k; */
  /* real<lower = -23, upper = -1.7> a; */
  vector<lower = -0.0001, upper = 0.0001>[n_ind] U_dos;
  vector<lower = -0.0001, upper = 0.0001>[n_ind] U_k;
  vector<lower = -0.0001, upper = 0.0001>[n_ind] U_a;
  real<lower=0, upper = 100> sigma;
}

model{
  
  //priors on random effects
  U_dos ~ normal(0, .0001);
  U_k ~ normal(0, .0001);
  U_a ~ normal(0, .0001);
  
  for(n in 1:n_obs){
	
	int idx = n_timepts[n];
	int n_meas = max(Time[n][:idx]);
	vector[n_meas] pred_fev1 = experimentFEV1(Cm[n][:idx],
										   Ve[n][:idx], Time[n][:idx],
										   exp(dos + U_dos[ind[n]]),
										   exp(k + U_k[ind[n]]),
										   -1 * exp(a + U_a[ind[n]]));
	
	int comp_idx[n_dFEV1[n]] = Time[n][dFEV1_measure_idx[n][:n_dFEV1[n]]];
	print("pred: ", pred_fev1[comp_idx] * -1);
	print("obs: ", dFEV1[n][:n_dFEV1[n]]);
	print("dos: ", exp(dos + U_dos[ind[n]]));
	print("k :", exp(k + U_k[ind[n]]));
	print("a :", exp(a + U_a[ind[n]]) * -1);
	// Likelihood
	dFEV1[n][:n_dFEV1[n]] ~ normal(pred_fev1[comp_idx] * -1, sigma);
	if(n == n_obs)
	  print("__lp: " , target());
  }
  
}

generated quantities{
  real log_lik = 0;
  real aic;
  for(n in 1:n_obs){
	int idx = n_timepts[n];
	vector[idx] pred_fev1 = experimentFEV1(Cm[n][:idx],
							   Ve[n][:idx], Time[n][:idx],
							   dos, k, a);

	int comp_idx[n_dFEV1[n]] = dFEV1_measure_idx[n][:n_dFEV1[n]];

	// Likelihood
	log_lik += normal_lpdf(dFEV1[n][:n_dFEV1[n]] | pred_fev1[comp_idx], sigma);
  }
  aic = 2 * 3 - 2 * log_lik;
}
