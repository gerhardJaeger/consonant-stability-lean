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

d <- read_tsv("../families/Indo-European/series_counts.tsv")

tree <- read.nexus("../families/Indo-European/phylogeny/heggartyetal2023.tree")

conversion <- read_csv("../families/Indo-European/phylogeny/IE_conversion.csv")

d <- inner_join(d, conversion, by = "glottocode")

# rename "series fullness" to "series_fullness"

d <- d %>% rename(series_fullness = `series fullness`) %>%
    select(name, series_fullness)

# prune tree to only include tips that are in d$name

tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d$name])

# check that all tips in d$name are in tree$tip.label
(sort(tree$tip.label) != sort(d$name)) %>% sum

# check for normality

d %>% ggplot(aes(x = series_fullness)) +
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
        edges <- rbind(edges, c(mother_id, daughter_id))
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
    if ((i %in% edges[,2]) && !(i %in% undefined)) {
        po <- c(po, i)
    }
}
edge_po <- match(po, edges[,2])

##

Y <- d$series_fullness

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

fit_pooling <- rstan::sampling(model_pooling, data = stan_data, chains = 4, iter = 2000, thin = 1)

print(fit_pooling, pars = c("mu", "sigma", "lp__"))

#       mean se_mean   sd  2.5%  25%  50%  75% 97.5% n_eff Rhat
# mu    3.50    0.00 0.07  3.37 3.46 3.50 3.55  3.64  3561    1
# sigma 0.59    0.00 0.05  0.50 0.56 0.59 0.63  0.71  3362    1
# lp__  1.71    0.02 1.02 -1.06 1.34 2.03 2.42  2.69  1798    1

(loo_pooling <- loo(fit_pooling))

#          Estimate   SE
# elpd_loo    -61.3  6.4
# p_loo         2.1  0.5
# looic       122.5 12.8
# ------
# Monte Carlo SE of elpd_loo is 0.0.
##

stan_data$Nnodes <- nrow(nodes)
stan_data$Nedges <- length(edge_po)
stan_data$observations <- match(d$name, nodes$label)
stan_data$root <- roots
stan_data$internal <- unique(edges[,1])
stan_data$Ninternal <- length(unique(edges[,1]))
stan_data$edges <- edges[edge_po,]
stan_data$edge_lengths <- edge_lengths[edge_po]


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

fit_ou <- rstan::sampling(model_ou, data = stan_data, chains = 4, iter = 10000, thin = 1)

saveRDS(fit_ou, "indo-european_ou.rds")

print(fit_ou, pars=c("mu", "sigma", "lambda", "lp__"))

#          mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# mu       3.55    0.00  0.12   3.31   3.47   3.54   3.62   3.80 16148    1
# sigma    0.65    0.00  0.13   0.49   0.57   0.63   0.71   0.93  1772    1
# lambda   0.62    0.01  0.44   0.24   0.44   0.56   0.73   1.30  1707    1
# lp__   -64.30    0.26 11.97 -90.82 -71.54 -63.30 -55.95 -43.85  2140    1


(loo_ou = loo(fit_ou))

#          Estimate   SE
# elpd_loo    -46.0  7.9
# p_loo        17.5  4.2
# looic        92.0 15.8
# ------

print(loo_compare(loo_pooling, loo_ou))

#        elpd_diff se_diff
# model2   0.0       0.0  
# model1 -15.2       7.1  