library(tidyverse)
library(gridExtra)
library(stringr)
library(reshape2)
library(rstan)
library(bridgesampling)
library(loo)
library(ape)
library(bayesplot)
library(coda)

options(mc.cores = parallel::detectCores())
rstan_options("auto_write" = TRUE)

d <- read_tsv("data/series_counts.tsv")

tree <- read.nexus("data/families/Pama-Nyungan/phylogeny/bouckaertetal2018.txt") %>% multi2di()

conversion <- read_csv("data/families/Pama-Nyungan/phylogeny/PN_conversion.csv")

d <- inner_join(d, conversion, by = "glottocode") %>%
    select(name, series_markedness_fullness)

# prune tree to only include tips that are in d$name

tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d$name])

# check that all tips in d$name are in tree$tip.label
(sort(tree$tip.label) != sort(d$name)) %>% sum

# check for normality

d %>% ggplot(aes(x = series_markedness_fullness)) +
    geom_density()


##


tips <- tree$tip.label

d <- d[match(tips, d$name), ]

##

node_matrix <- c()
node_number <- 1
tn <- 1

tr_nodes <- length(tree$tip.label) + tree$Nnode
for (i in 1:tr_nodes) {
  tl <- ""
  if (i <= length(tree$tip.label)) {
    tl <- tree$tip.label[i]
  }
  node_matrix <- rbind(node_matrix, c(node_number, tn, i, tl))
  node_number <- node_number+1
}


nodes <- tibble(
    node_number = as.numeric(node_matrix[,1]),
    tree_number = as.numeric(node_matrix[,2]),
    internal_number = as.numeric(node_matrix[,3]),
    label = node_matrix[,4]
)

get_global_id <- function(tree_number, local_id) which((nodes$tree_number == tree_number) & (nodes$internal_number == local_id))


edges <- c()
edge_lengths <- c()
tn <- 1
tr <- tree
if (length(tr$edge) > 0 ) {
    for (en in 1:nrow(tr$edge)) {
        mother <- tr$edge[en, 1]
        daughter <- tr$edge[en, 2]
        mother_id <- get_global_id(tn, mother)
        daughter_id <- get_global_id(tn, daughter)
        el <- tr$edge.length[en]
        newrow <- c(mother_id, daughter_id)
        edges <- rbind(edges, newrow)
        edge_lengths <- c(edge_lengths, el)
    }
}


roots <- sort(setdiff(edges[,1], edges[,2]))
tips <- which(nodes$label != "")


get_postorder <- function(roots, edges) {
    input <- roots
    output <- c()
    while (length(input) > 0) {
        pivot <- input[1]
        daughters <- setdiff(edges[edges[,1] == pivot, 2], output)
        if (length(daughters) == 0) {
            input <- input[-1]
            output <- c(output, pivot)
        } else {
            input <- c(daughters, input)
        }
    }
    return(output)
}


undefined <- c()
po <- c()
for (i in get_postorder(roots, edges)) {
    if ((i %in% edges[,2]) & !(i %in% undefined)) {
        po <- c(po, i)
    }
}
edge_po <- match(po, edges[,2])

##

Y <- d$series_markedness_fullness

stan_data <- list()
stan_data$Y <- Y
stan_data$Ndata <- length(Y)

pooling_code <- "
data {
    int<lower=1> Ndata;
    vector[Ndata] Y;
}

parameters {
    real mu;
    real<lower=0> sigma;
}

model {
    Y ~ normal(mu, sigma);
    mu ~ normal(0, 10);
    sigma ~ cauchy(0, 5);
}

generated quantities {
    vector[Ndata] log_lik;
    for (i in 1:Ndata) {
        log_lik[i] = normal_lpdf(Y[i] | mu, sigma);
    }
}
"


model_pooling <- stan_model(model_code = pooling_code)

fit_pooling <- rstan::sampling(
  model_pooling,
  data = stan_data,
  chains = 4,
  iter = 2000,
  thin = 1
)

print(fit_pooling, pars = c("mu", "sigma", "lp__"))

#        mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# mu     1.21    0.00 0.03  1.15  1.19  1.21  1.23  1.27  3661    1
# sigma  0.45    0.00 0.02  0.41  0.44  0.45  0.47  0.50  3710    1
# lp__  63.39    0.02 0.98 60.76 62.98 63.68 64.11 64.37  1997    1

(loo_pooling <- loo(fit_pooling))
#          Estimate   SE
# elpd_loo   -135.9 15.2
# p_loo         3.1  0.6
# looic       271.7 30.5

##

stan_data$Nnodes <- nrow(nodes)
stan_data$Nedges <- length(edge_po)
stan_data$observations <- match(d$name, nodes$label)
stan_data$root <- roots
stan_data$internal <- unique(edges[, 1])
stan_data$Ninternal <- length(unique(edges[, 1]))
stan_data$edges <- edges[edge_po, ]
stan_data$edge_lengths <- edge_lengths[edge_po] + 1e-6



ou_code <- "
data {
    int<lower=1> Ndata;
    vector[Ndata] Y;
    int<lower=1> Nnodes;
    int<lower=1> Nedges;
    int<lower=1> Ninternal;
    int<lower=1,upper=Nnodes> root;
    int<lower=1,upper=Nnodes> observations[Ndata];
    int<lower=1,upper=Nnodes> internal[Ninternal];
    int<lower=1, upper=Nnodes> edges[Nedges, 2];
    vector[Nedges] edge_lengths;
}

parameters {
    vector[Ninternal] z_internal;
    real mu;
    real<lower=0> sigma;
    real<lower=0> lambda;
    }

transformed parameters {
   vector[Nnodes] z = rep_vector(0, Nnodes);
   for (i in 1:Ndata) {
    z[observations[i]] = Y[i];
   }
   for (i in 1:Ninternal) {
    z[internal[i]] = z_internal[i];
   }
}

model {
    mu ~ normal(0, 10);
    sigma ~ exponential(1);
    lambda ~ cauchy(0, 2.5);
    
    target += normal_lpdf(z[root] | mu, sigma * sqrt(2 * lambda));
    
    for (i in 1:Nedges) {
        int j = Nedges - i + 1;
        int mother_node = edges[j, 1];
        int daughter_node = edges[j, 2];
        real local_mu = mu + (z[mother_node] - mu) * exp(-lambda * edge_lengths[j]);
        real local_sigma = sigma * sqrt((1 - exp(-2 * lambda * edge_lengths[j])) / (2 * lambda));
        target += normal_lpdf(z[daughter_node] | local_mu, local_sigma);
    }
}

generated quantities {
   vector[Nnodes] log_lik_z;
   vector[Ndata] log_lik;

   log_lik_z[root] = normal_lpdf(z[root] | mu, sigma * sqrt(2 * lambda));
   for (i in 1:Nedges) {
        int j = Nedges - i + 1;
        int mother_node = edges[j, 1];
        int daughter_node = edges[j, 2];
        real local_mu = mu + (z[mother_node] - mu) * exp(-lambda * edge_lengths[j]);
        real local_sigma = sigma * sqrt((1 - exp(-2 * lambda * edge_lengths[j])) / (2 * lambda));
        log_lik_z[daughter_node] = normal_lpdf(z[daughter_node] | local_mu, local_sigma);
    }
    for (i in 1:Ndata) {
        log_lik[i] = log_lik_z[observations[i]];
    }
}

"


model_ou <- stan_model(model_code = ou_code)


fit_ou <- rstan::sampling(
  model_ou,
  data = stan_data,
  chains = 4,
  iter = 10000,
  thin = 10,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

saveRDS(fit_ou, "pama-nyungan_ou.rds")

print(fit_ou, pars = c("mu", "sigma", "lambda", "lp__"))


#           mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
# mu        1.21    0.00  0.03    1.16    1.19    1.22    1.23    1.27  1991    1
# sigma     1.09    0.02  0.65    0.28    0.66    0.96    1.34    2.65  1584    1
# lambda    4.01    0.20  6.81    0.20    1.11    2.26    4.40   17.80  1137    1
# lp__   -270.25    0.35 14.62 -299.35 -279.98 -269.94 -259.98 -243.48  1713    1

(loo_ou <- loo(fit_ou))

#          Estimate   SE
# elpd_loo   -135.8 15.4
# p_loo         3.0  0.6
# looic       271.5 30.7


print(loo_compare(loo_pooling, loo_ou))


#        elpd_diff se_diff
# model2  0.0       0.0   
# model1 -0.1       0.