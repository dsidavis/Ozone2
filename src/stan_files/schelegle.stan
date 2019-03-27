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
  
  vector UOS(real DR, vector CumFrDos, real Dos, vector t) {
	int n = dims(t)[1];
	vector[n] UOS;
	for(i in 1:n)
	  UOS[i] = DR / (1 + exp(-20 * (t[i] - (t[i] / CumFrDos[i]))));
	
	return UOS;
  }
  
  vector deltaX(vector uos, real k, real x_prev){
	int n = dims(uos)[1];
	real r = (1 - exp(-k));
	vector[n] ceoss = uos / k;
	vector[n] x = rep_vector(0.0, n);

	x[1] = x_prev;

	if(n > 1)
	  for(i in 2:n)
		x[i] = x[i-1] + ((ceoss[i-1] - x[i-1]) * r);
	
	return x;
  }

  vector experimentFEV1(vector o3, vector ve, int[] t_stop,
						real dos, real k, real a){
	int m = dims(o3)[1];
	int n = max(t_stop);
	vector[m] dFEV1 = rep_vector(0, m);
	real x_previous = 0;
	real dr = 0;
	real FrDos;
	real FrDos_previous = 0;
	
	for(i in 1:(m - 1)){
	  int n_t = t_stop[i+1] - t_stop[i]+1;
	  vector[n_t] t = seq(t_stop[i], t_stop[i+1]);
	  vector[n_t] CumFrDos;
	  vector[n_t] uos;
	  vector[n_t] x;
	  dr = o3[i] * ve[i] * 1.96;
	  FrDos = dr / dos;
	  for(j in 1:n_t)
		CumFrDos[j] = FrDos_previous + FrDos * j;
	  /* print(CumFrDos); */
	  uos = UOS(dr, CumFrDos, dos, t);
	  /* print("UOS: ", uos); */
	  x = deltaX(uos, k, x_previous);
	  FrDos_previous = CumFrDos[n_t];
	  /* print("x = ", x); */
	  x_previous = x[n_t];
	  dFEV1[i] = x_previous;		
	}
	/* print("dFEV1 = ", dFEV1); */
	return dFEV1 * a;
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
  real<lower = 0> sigma;
}
transformed data{
  /* real k = 0.02; */
  /* real dos = 900; */
  /* real a = -0.02; */
}
parameters{
  /* real dos; */
  /* real k; */
  /* real a; */
  /* real<lower=0> sigma; */
}

model{

  /* dos ~ normal(1000, 1000); */
  /* k ~ normal(0,1); */
  /* a ~ normal(0,1); */
  
  for(n in 1:n_obs){
	int idx = n_timepts[n];
	int n_meas = max(Time[n][:idx]);
	vector[idx] pred_fev1 = experimentFEV1(Cm[n][:idx],
							   Ve[n][:idx], Time[n][:idx],
							   dos, k, a);
	
	int comp_idx[n_dFEV1[n]] = dFEV1_measure_idx[n][:n_dFEV1[n]];
	/* print("pred: ", pred_fev1[comp_idx] * -100); */
	/* print("obs: ", dFEV1[n][:n_dFEV1[n]]); */
	// Likelihood
	dFEV1[n][:n_dFEV1[n]] ~ normal(pred_fev1[comp_idx] * -100, sigma);
	/* print("__lp: " , target()); */
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
