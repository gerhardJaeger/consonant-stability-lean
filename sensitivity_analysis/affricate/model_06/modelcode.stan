data {
  int<lower=1> Ndata;
  int<lower=1> Nnodes;
  int<lower=1> Nedges;
  vector[Ndata] y; 
  array[Nedges, 2] int<lower=1, upper=Nnodes> edges;      // parent → child
  vector<lower=0>[Nedges] edge_lengths;                   // edge lengths
  int<lower=1, upper=Nnodes> root_node;                   // index of root node
}

transformed data {
   vector[Ndata] y_transformed;
   for (i in 1:Ndata) {
    // y_transformed[i] = log(y[i] + 0.001);
     y_transformed[i] = y[i];
   }
}

parameters {
  vector[Nnodes] z_std;       // standard normal reparam for latent z
//   real<lower=0> sigma;        // OU diffusion
  real<lower=0> V;
  real<lower=0> t_half;       // phylogenetic half-life
  real mu;           // OU mean
  real<lower=0> tau;          // measurement noise
}

transformed parameters {
  vector[Nnodes] z;
  real lambda = log(2) / t_half; // ensure lambda is positive
  real sigma = sqrt(2 * lambda * V);

  // root node
  z[root_node] = mu + (sigma / sqrt(2 * lambda)) * z_std[root_node];

  // recursive evolution
  for (e in 1:Nedges) {
    int edge_index = Nedges - e + 1; // reverse order for recursion
    int parent = edges[edge_index, 1];
    int child = edges[edge_index, 2];
    real len = edge_lengths[edge_index];

    real decay = exp(-lambda * len);
    real s = sigma * sqrt(-expm1(-2 * lambda * len) / (2 * lambda));
    real mn = mu + (z[parent] - mu) * decay;

    z[child] = mn + s * z_std[child];
  }
}

model {
  // Priors
  // sigma ~ lognormal(0, 1);
  V ~ lognormal(0, 1);
  // t_half ~ lognormal(0, 1);
  t_half ~ lognormal(2, 0.5);
  // mu ~ normal(0, 2);
  mu ~ normal(2, 1);
  //tau ~ exponential(10);
  tau ~ lognormal(-1, 0.5);

  z_std ~ normal(0, 1);  // standard normal prior

  // Likelihood
  for (i in 1:Ndata) {
    y_transformed[i] ~ normal(z[i], tau); // measurement noise
  }
}

generated quantities {
  vector[Ndata] log_lik;
  vector[Ndata] y_rep;     // posterior predictive
  vector[Ndata] y_prior;   // prior predictive
  real eps = 0.001;
  // Posterior predictive samples (uses fitted alpha, beta, z)
  for (i in 1:Ndata) {
    log_lik[i] = normal_lpdf(y_transformed[i] | z[i], tau);
    // y_rep[i] = exp(normal_rng(z[i], tau))-0.001;
    y_rep[i] = normal_rng(z[i], tau);
  }

  // Prior predictive samples: draw z_prior from prior
  {
    
    vector[Nnodes] z_prior_std;
    for (j in 1:Nnodes)
      z_prior_std[j] = normal_rng(0, 1);

    // real mu_prior = normal_rng(0, 2);
    real mu_prior = normal_rng(2, 1);
    // real t_half_prior = lognormal_rng(0, 1);
    real t_half_prior = lognormal_rng(2, 0.5);
    // real sigma_prior = lognormal_rng(0, 1);
    real V_prior = lognormal_rng(0, 1);
    real lambda_prior = log(2) / t_half_prior;
    real sigma_prior = sqrt(2 * lambda_prior * V_prior);
    // real tau_prior = exponential_rng(10);
    real tau_prior = lognormal_rng(-1, 0.5);
    vector[Nnodes] z_prior;
    z_prior[root_node] = mu_prior + (sigma_prior / sqrt(2 * lambda_prior)) * z_prior_std[root_node];

    for (e in 1:Nedges) {
      int edge_index = Nedges - e + 1;
      int parent = edges[edge_index, 1];
      int child = edges[edge_index, 2];
      real len = edge_lengths[edge_index];
      real decay = exp(-lambda_prior * len);
      real s = sigma_prior * sqrt(-expm1(-2 * lambda_prior * len) / (2 * lambda_prior));
      real mn = mu_prior + (z_prior[parent] - mu_prior) * decay;
      z_prior[child] = mn + s * z_prior_std[child];
    }

    for (i in 1:Ndata) {
      //y_prior[i] = exp(normal_rng(z_prior[i], tau_prior))-0.001;
      y_prior[i] = normal_rng(z_prior[i], tau_prior);
    }
  }
}
