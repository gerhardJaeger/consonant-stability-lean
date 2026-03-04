## Uralic family OU analysis

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
library(codetools)
library(ggtree)

options(mc.cores = parallel::detectCores())
rstan_options("auto_write" = TRUE)

fam <- "../families/Uralic/"
d <- read_tsv(str_c(fam, "series_counts.tsv"))

tree <- read.nexus(str_c(fam, "phylogeny/uralic-samoyedic_constraint.mcc.nex")) %>% 
    multi2di()

conversion <- read_csv(str_c(fam, "phylogeny/uralic_conversion.csv"))

d <- inner_join(d, conversion, by = "glottocode")

# rename "series fullness" to "series_fullness"

d <- d %>% rename(series_markedness_fullness = `series markedness fullness`) %>%
    select(name, series_markedness_fullness)

# prune tree to only include tips that are in d$name

tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d$name])

# check that all tips in d$name are in tree$tip.label
(sort(tree$tip.label) != sort(d$name)) %>% sum

d %>% ggplot(aes(x = series_markedness_fullness)) +
    geom_density()

tips <- tree$tip.label

d <- d[match(tips, d$name), ]

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

(loo_pooling <- loo(fit_pooling))

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
    
    target += normal_lpdf(z[root] | mu, sigma / sqrt(2 * lambda));
    
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

   log_lik_z[root] = normal_lpdf(z[root] | mu, sigma / sqrt(2 * lambda));
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
saveRDS(fit_ou, "uralic_ou.rds")

print(fit_ou, pars = c("mu", "sigma", "lambda", "lp__"))

(loo_ou <- loo(fit_ou))

print(loo_compare(loo_pooling, loo_ou))

value_df <- tibble(series=apply(extract(fit_ou, pars="z")$z, 2, mean)) %>%
    mutate(node=1:nrow(.))

merged_df <- fortify(tree) %>%
    merge(value_df, by="node")

options(repr.plot.width=20, repr.plot.height=15)
ggtree(tree) +
  geom_tiplab() +  # Add tip labels
  geom_point(data = merged_df, aes(color = series), size = 5) +  # Plot continuous values
  scale_color_continuous(low = "green", high = "red") +  # Define color scale
  geom_text(data = merged_df, aes(label = round(series, 2))) +  # Add text labels for series at nodes
  theme_tree2()

