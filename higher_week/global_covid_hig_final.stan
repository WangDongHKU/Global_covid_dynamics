functions {
    vector fmax_vec(vector a, vector b) {
     return 0.5 * (a + b + fabs(a - b));
  }
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
      //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
    b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}
data {
  int<lower=0> K; // number of PROVIENCES
  int<lower=0> W; // number of weeks
  vector [K] N; // population size
  int<lower=0> cases[K, W];    // ILI+ PROXY    //
  matrix [K,1724] NPI1;    
  matrix[K,K] Mobility[1709];
  real error;
  int num_data;             // number of data points
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  real X[num_data];
  row_vector[num_knots+spline_degree-1] a_raw1[K];
  #matrix [K,W] drift;
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_data] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  vector[K] zerovector = rep_vector(0,K);
  real gamma = 1.0/10.0;
  real q = 1.0/610;
  real sigma = 1.0/5;
ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis){
  B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));}
B[num_knots + spline_degree - 1, num_data] = 1;
}

parameters {
  vector <lower=0,upper=1>[K] thetap;
  real<lower=0.96,upper=1> S0[2];
  real<lower=0,upper=0.3> E0[2];
  real<lower=0,upper=0.3> I0[2];
  real <lower=0,upper=1> npi;
  real <lower=0>phi[K];
  row_vector<upper=10>[num_basis] a_raw[K];
 real <lower=0> rate;
 vector <lower=0,upper=1>[K] p;
 vector <lower=0> [3] v;
 vector <lower=1> [K] sp;
 //real drift_rate;
}

transformed parameters {
  matrix <lower=0>[K,W*7]S;
  matrix <lower=0>[K,W*7]E;
  matrix <lower=0>[K,W*7]I;
  matrix <lower=0>[K,W*7]R;
  matrix[K,W] PILI;  // weekly simulated data

  matrix<lower=0>[K,W*7]Rt1;

  vector[K]newS;
  vector[K]newE;
  vector[K]newI;
  vector[K]newR;
  matrix[K,num_data] beta1;

  vector[K] total;
  vector[K] npifactor;
  vector[K] move_vec;
  vector[K] moveE_vec;
  vector[K] ptemp;
  real vtemp;
  //vector[K] sptemp;

  row_vector[num_basis] a[K];
  matrix[K,num_data] beta1_KT;

   PILI = rep_matrix(0, K, W);
  for (k in 1:K){
  a[k,] = a_raw[k,];
  beta1_KT[k,] =  exp(to_row_vector(a[k,]*B));
  }
  //print("beta1_KT[k,] :", beta1_KT[,15]);
  // Initial
  for (k in 1:K) {
    beta1[k,1]=0.5;
    Rt1[k,1]= 3;
    if(k==1){
    S[k,1] = S0[1]; 
    E[k,1] = E0[1];
    I[k,1] = I0[1];
    R[k,1] = 0;}
    else{
    S[k,1] = S0[2]; 
    E[k,1] = E0[2];
    I[k,1] = I0[2];
    R[k,1] = 0;
    }
    total[k] = S[k,1]+E[k,1]+I[k,1]+R[k,1]+ 1e-10;
    S[k,1] = S[k,1]/total[k];
    E[k,1] = E[k,1]/total[k];
    I[k,1] = I[k,1]/total[k];
    R[k,1] = R[k,1]/total[k];
  }
    PILI[,1] = sigma*E[,1];  
  
  // SIR
for (w in 1:W) {
    vtemp = v[1] + (w >= 52) * (v[2] - v[1]) + (w >= 104) * (v[3] - v[2]);
    //sptemp = sp[,1] + (w >= 52) * (sp[,2] - sp[,1]) + (w >= 104) * (sp[,3] - sp[,2]);
    ptemp = p;//[1] + (w >= 52) * (p[2] - p[1]) + (w >= 104) * (p[3] - p[2]);

  for (t in max(2,(w-1)*7+1):(w*7)) {
    // 计算传播率
    beta1[,t] = beta1_KT[:,t] .* exp(-npi* NPI1[:,t]);
    //{print("beta1[,t] :", beta1[,t]);}
    // 计算NPI影响因子，添加数值稳定性保护
    npifactor = 1 - thetap ./ (1 + exp(-rate * NPI1[:,t])); // NPI factor vector

    move_vec = ((Mobility[t-1])* I[:, t-1]  - (rep_row_vector(1, K)*(Mobility[t-1]))' .* I[:, t-1]) .* npifactor ./ N;
    moveE_vec = ((Mobility[t-1])* E[:, t-1] - (rep_row_vector(1, K)*(Mobility[t-1]))'.*E[:, t-1]) ./ N ;

    // 修改状态更新计算，添加数值稳定性保护
    newS = (q * R[:, t-1] - (1-ptemp) .* S[:, t-1] .* beta1[, t] .* (I[:, t-1] + move_vec) - ptemp .* sp .* S[:, t-1] .* beta1[, t] .* (I[:, t-1] + move_vec  + vtemp .* square(I[:, t-1] + move_vec ))); 
    S[:, t] = fmax_vec(zerovector, S[:, t-1] + newS);
 
    newE = ((1-ptemp) .* S[:, t-1] .* beta1[, t] .* (I[:, t-1] + move_vec ) + ptemp .* sp .* S[:, t-1] .* beta1[, t] .* (I[:, t-1] + move_vec  + vtemp .* square(I[:, t-1] + move_vec )) + moveE_vec  - sigma * E[:, t-1]); 
    E[:, t] = fmax_vec(zerovector, E[:, t-1] +  newE);

    newI = (sigma * E[:, t-1] - gamma * I[:, t-1]);
    I[:, t] = fmax_vec(zerovector, I[:, t-1] +  newI);

    newR = (gamma * I[:, t-1] - q * R[:, t-1]);
    R[:, t] = fmax_vec(zerovector, R[:, t-1] +  newR);

    // 修改归一化过程，添加数值稳定性保护
    total = S[:, t] + E[:, t] + I[:, t] + R[:, t]+ 1e-10;

    S[:, t] = S[:, t] ./ total;
    E[:, t] = E[:, t] ./ total;
    R[:, t] = R[:, t] ./ total;
    I[:, t] = I[:, t] ./ total;

    PILI[:, w] = PILI[:, w]  +  sigma * E[:, t] .* N .* (1-npifactor);
    Rt1[:, t] = S[:, t] .*  beta1[, t] / gamma;
}
}
}

model {
  // 参数先验
  npi ~ normal(0,0.001);   
  S0[1] ~ lognormal(log(0.982),0.01);
  E0[1] ~ lognormal(log(0.0005),0.001);
  I0[1] ~ lognormal(log(0.0005),0.001);
  S0[2] ~ lognormal(log(0.992),0.01);
  E0[2] ~ lognormal(log(0.0001),0.001);
  I0[2] ~ lognormal(log(0.0003),0.001);
  rate ~ normal(0,0.05);

    v ~ normal(0,1);
    p ~ normal(0.0,0.5);
    sp ~ normal(1,2);
    phi ~ normal(0,3);
    thetap ~ normal(0.02591,0.1);
    //drift_rate ~ normal(0,0.05);
  for (k in 1:K){
    a_raw[k,] ~ normal(a_raw1[k,], 0.001); 
    cases[k,]  ~  neg_binomial_2(PILI[k,] + error,  phi[k]); 
    }
}

generated quantities {
  real pred_cases[K,W];
  real log_lik[K,W];

for (k in 1:K){ 
    for (w in 1:W) { 
	  pred_cases[k,w] = neg_binomial_2_rng(PILI[k,w] + error, phi[k]);
	  log_lik[k,w] = neg_binomial_2_lpmf(cases[k,w] | PILI[k,w] + error,  phi[k]);
}
}
} 
