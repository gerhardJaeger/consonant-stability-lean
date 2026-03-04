library(BBMV)
library(geiger)
library(tidyverse)
library(coda)
library(httr)


tree <- read.nexus("global-language-tree-MCC-labelled.tree")




taxa <- tree$tip.label
# split each taxon at the first underscore and keep the first part
taxa <- sapply(strsplit(taxa, "_"), function(x) x[1])
tree$tip.label <- taxa

d <- read_tsv("../charts/series_counts.tsv") %>%
    filter(glottocode %in% taxa)

## get the Glottolog data and add Family_ID to d

url <- "https://zenodo.org/record/10804582/files/glottolog/glottolog-cldf-v5.0.zip?download=1"
zip_path <- "glottolog-cldf-v5.0.zip"

GET(url, write_disk(zip_path, overwrite = TRUE))

unzip(zip_path, exdir = ".")
unzipped_files <- list.files("glottolog-glottolog-cldf-4dbf078", recursive = TRUE)
languages_file <- unzipped_files[grep("languages.csv", unzipped_files)]

languages <- read_csv(file.path("glottolog-glottolog-cldf-4dbf078", languages_file))

file.remove(zip_path)
unlink("glottolog-glottolog-cldf-4dbf078", recursive = TRUE)



family_ids <- languages %>%
    filter(Glottocode %in% d$glottocode) %>%
    drop_na(Family_ID) %>%
    group_by(Family_ID) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    select(Family_ID) %>%
    pull

families <- languages$Name[match(family_ids, languages$ID)]
names(families) <- family_ids

d <- d %>%
    inner_join(
        select(languages, Glottocode, Family_ID) %>% drop_na %>% distinct %>% rename(glottocode = Glottocode),
        by = "glottocode"
    ) %>%
    mutate(Family = as.character(families[Family_ID]))
##

# prune tree to d$glottocode

tree <- drop.tip(tree, setdiff(tree$tip.label, d$glottocode))

d <- d[match(tree$tip.label, d$glottocode), ]
##
l_trait <- log(d$series_markedness_fullness)
names(l_trait) <- d$glottocode

d <- d %>%
    mutate(log_series_markedness_fullness = log(series_markedness_fullness))
##

tsv_string <- capture.output(write_tsv(d[,-c(14, 15)], stdout()))
tsv_string <- paste(tsv_string[-1], collapse = "\n")

nex <- paste(
    "#NEXUS\nBegin data;\nDimensions ntax=",
    length(tree$tip.label),
    " nchar=",
    ncol(d) - 3,
    ";\nFormat datatype=Continuous missing=? gap=-;\nMatrix\n",
    tsv_string,
    "\n;\nEnd;\n",
    sep = "",
    collapse = ""
)

cat(nex, file = "../data/data_pruned.nex")


write_csv(d, "../data/data_pruned.csv")

write.nexus(tree, file = "../data/tree_pruned.nex", translate = TRUE)


##

ll_FPK4  <- lnL_FPK(tree, l_trait, Npts = 50, a = NULL, b = NULL, c = NULL) # the full model

ll_FPK2 <- lnL_FPK(tree, l_trait, Npts = 50, a = 0, b = NULL, c = NULL)

ll_FPK0 <- lnL_FPK(tree, l_trait, Npts = 50, a = 0, b = 0, c = 0)

fit4 <- find.mle_FPK(model = ll_FPK4)
get.landscape.FPK(fit = fit4)


fit2 <- find.mle_FPK(model = ll_FPK2)
get.landscape.FPK(fit = fit2) # this shape of the landscape cannot have 2 peaks


fit0 <- find.mle_FPK(model = ll_FPK0)
get.landscape.FPK(fit = fit0) # this one is forced to be flat


fit4$aic
fit2$aic
fit0$aic

charac_time(fit = fit4)
charac_time(fit = fit2)
charac_time(fit = fit0)

max(branching.times(tree))

summary(tree$edge.length)

plot(fit4$root, type = "l")

fit4$par

Uncertainty_FPK(
    fit = fit4,
    tree,
    trait = l_trait,
    Npts = 50,
    effort_uncertainty = 100,
    scope_a = c(-5, 10),
    scope_b = c(-5, 5),
    scope_c = c(-2, 2)
)

OU <- fitContinuous(phy = tree, dat = l_trait, model = "OU")
OU$opt$lnL
fit2$lnL

BM <- fitContinuous(phy = tree, dat = l_trait, model = "BM")
BM$opt$aic

ACE_nodes <- ACE_FPK(fit4, specific.point = NULL)

plot(ACE_nodes[[length(tree$tip.label) + 1]], type = "l")

bounds <- log(c(1, 25))

ll_BBMV4 <- lnL_BBMV(tree, l_trait, Npts = 50, bounds = bounds, a = NULL, b = NULL, c = NULL)
ll_BBMV2 <- lnL_BBMV(tree, l_trait, Npts = 50, bounds = bounds, a = 0, b = NULL, c = NULL)
ll_BBMV1 <- lnL_BBMV(tree, l_trait, Npts = 50, bounds = bounds, a = 0, b = 0, c = NULL)
ll_BBMV0 <- lnL_BBMV(tree, l_trait, Npts = 50, bounds = bounds, a = 0, b = 0, c = 0) # this is the BBM model

fit4b <- find.mle_FPK(model = ll_BBMV4)
get.landscape.FPK(fit = fit4b)

fit2b <- find.mle_FPK(model = ll_BBMV2)
get.landscape.FPK(fit = fit2b)

fit1b <- find.mle_FPK(model = ll_BBMV1)
get.landscape.FPK(fit = fit1b)

fit0b <- find.mle_FPK(model = ll_BBMV0)
get.landscape.FPK(fit = fit0b)

fit4b$aic
fit2b$aic
fit1b$aic
fit0b$aic


Uncertainty_FPK(
    fit = fit4b,
    tree,
    trait = l_trait,
    Npts = 50,
    effort_uncertainty = 100,
    scope_a = c(-10, 10),
    scope_b = c(-10, 10),
    scope_c = c(-10, 10)
)

par(mfrow = c(1, 1))
ACE_nodes <- ACE_FPK(fit4b, specific.point = NULL)
plot(ACE_nodes[[length(tree$tip.label) + 1]], type = "l")

par(mfrow = c(1, 1))
plot(fit4b$root, type = "l")

charac_time(fit = fit4b)

##

fm2 <- d %>%
    select(glottocode, Family) %>%
    group_by(Family) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    select(Family) %>%
    pull

trees <- list()
TRAITS <- list()


for (fm in fm2) {
    d_fm <- d %>%
        filter(Family == fm) %>%
        select(glottocode, series_markedness_fullness) %>%
        drop_na %>%
        mutate(series_markedness_fullness = log(series_markedness_fullness))
    fm_trait <- d_fm$series_markedness_fullness
    names(fm_trait) <- d_fm$glottocode
    TRAITS[[fm]] <- fm_trait
    trees[[fm]] <- drop.tip(tree, setdiff(tree$tip.label, d_fm$glottocode))
}
##

lnls <- ks <- rep(NA, length(trees))
fits <- list()
good_families <- c()
for (i in seq_along(trees)) {
    print(names(trees)[i])
    lnl_temp <- lnL_BBMV(trees[[i]], TRAITS[[i]], bounds = bounds, a = NULL, b = NULL, c = NULL, Npts = 50)
    fit_temp <- try(find.mle_FPK(model = lnl_temp, method = "Nelder-Mead", init.optim = NULL), silent = FALSE)
    if (!inherits(fit_temp, "try-error")) {
        fits[[i]] <- fit_temp
        names(fits)[i] <- paste("fit_clade_", i, sep = "")
        lnls[i] <- fit_temp$lnL
        ks[i] <- fit_temp$k
        good_families <- c(good_families, names(trees)[i])
    }
}
good_indices <- which(!is.na(lnls))

fitmFPK4 <- list(
    lnL = sum(lnls[good_indices]),
    aic = 2 * (sum(ks[good_indices]) - sum(lnls[good_indices])),
    k = sum(ks[good_indices]),
    fits = fits[good_indices]
)

lnl_FPK_multiclades_same_V_same_sig2 <- function(
    trees, traits, bounds = bounds, a = NULL, b = NULL, c = NULL, Npts = 50
) {
    if (length(trees) != length(traits)) {
        stop("The list of trees and the list of traits differ in length.")
    }
    if (length(trees) == 1) {
        stop("There is only one tree and trait vector: use the function lnl_BBMV instead")
    }
    for (i in 1:length(trees)) {
        if (sum(trees[[i]]$tip.label %in% names(traits[[i]])) <
            max(length(traits[[i]]), length(trees[[i]]$tip.label))) {
            stop(paste("Tip names in tree ", i, " do not match names of corresponding trait vector"))
        }
    }
    SEQ <- seq(from = -1.5, to = 1.5, length.out = Npts)
    trees_formatted <- list()
    for (i in 1:length(trees)) {
        trees_formatted[[i]] <- FormatTree_bounds(trees[[i]], traits[[i]], V = rep(0, Npts), bounds = bounds)
    }
    ncoeff <- (is.null(a) == T) + (is.null(b) == T) + (is.null(c) == TRUE) 
    if (is.null(a) == FALSE) {
        if (is.null(b) == F) {
            # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work
            if (is.null(c) == F) {
                fun_text <- "fun=function(X){return("
                for (i in 1:length(trees)) {
                    fun_text <- paste(
                        fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)",
                        sep = ""
                    )
                }
                fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
                fun <- eval(parse(text = fun_text))
            } else {
                fun_text <- "fun=function(X){return("
                for (i in 1:length(trees)) {
                    fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[2]*SEQ),bounds=bounds)", sep = "")
                }
                fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
                fun <- eval(parse(text = fun_text))
            }
        } else {
            fun_text <- "fun=function(X){return("
            for (i in 1:length(trees)) {
                fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+X[2]*SEQ^2+X[3]*SEQ),bounds=bounds)", sep = "")
            }
            fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
            fun <- eval(parse(text = fun_text))
        }
    } else {
        fun_text <- "fun=function(X){return("
        for (i in 1:length(trees)) {
            fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[1],dMat=DiffMat_backwards(X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ),bounds=bounds)", sep = "")
        }
        fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
        fun <- eval(parse(text = fun_text))
    }
    return(
        list(
            fun = fun,
            ncoeff = ncoeff,
            par_fixed = list(a = a, b = b, c = c, bounds = bounds),
            trees = trees,
            traits = traits,
            Npts = Npts
        )
    )
}

testFPK4 <- lnl_FPK_multiclades_same_V_same_sig2(
    trees = trees[good_indices],
    traits = TRAITS[good_indices],
    bounds = bounds,
    a = NULL, b = NULL, c = NULL, Npts = 50
)

fitFPK4 <- find.mle_FPK_multiple_clades_same_V_same_sig2(model = testFPK4, method = "Nelder-Mead", init.optim = NULL)



lnl_FPK_multiclades_same_V_different_sig2 <- function(trees, traits, bounds=bounds, a = NULL, b = NULL, c = NULL, Npts = 50) {
    if (length(trees) != length(traits)) {
        stop("The list of trees and the list of traits differ in length.")
    }
    if (length(trees) == 1) {
        stop("There is only one tree and trait vector: use the function lnl_FPK instead")
    }
    for (i in 1:length(trees)) {
        if (sum(trees[[i]]$tip.label %in% names(traits[[i]])) < max(length(traits[[i]]), length(trees[[i]]$tip.label))) {
            stop(paste("Tip names in tree ", i, " do not match names of corresponding trait vector"))
        }
    }
    SEQ <- seq(from = -1.5, to = 1.5, length.out = Npts)
    trees_formatted <- list()
    for (i in 1:length(trees)) {
        trees_formatted[[i]] <- FormatTree_bounds(trees[[i]], traits[[i]], V = rep(0, Npts), bounds = bounds)
    }
    ncoeff <- (is.null(a) == T) + (is.null(b) == T) + (is.null(c) == T) # OK
    npar <- ncoeff + length(trees)
    par_names <- c()
    for (i in 1:length(trees)) {
        par_names <- c(par_names, paste("dCoeff_tree_", i, sep = ""))
    }
    par_names <- c(par_names, "a", "b", "c")
    if (is.null(a) == F) {
        if (is.null(b) == F) {
            # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work
            if (is.null(c) == F) {
                fun_text <- "fun=function(X){return("
                for (i in 1:length(trees)) {
                    fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[", i, "],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)", sep = "")
                }
                fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
                fun <- eval(parse(text = fun_text))
            }

            # only c varies (e.g. flat landscape if a=b=0): seems to work
            else {
                fun_text <- "fun=function(X){return("
                for (i in 1:length(trees)) {
                    fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[", i, "],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[", length(trees) + 1, "]*SEQ),bounds=bounds)", sep = "")
                }
                fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
                fun <- eval(parse(text = fun_text))
            }
        }
        # only a is fixed (e.g. quadratic landscape if a=0): seems to work
        else {
            fun_text <- "fun=function(X){return("
            for (i in 1:length(trees)) {
                fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[", i, "],dMat=DiffMat_backwards(a*SEQ^4+X[", length(trees) + 1, "]*SEQ^2+X[", length(trees) + 2, "]*SEQ),bounds=bounds)", sep = "")
            }
            fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
            fun <- eval(parse(text = fun_text))
        }
    }

    # the full model: no parameter fixed
    else {
        fun_text <- "fun=function(X){return("
        for (i in 1:length(trees)) {
            fun_text <- paste(fun_text, "-LogLik_bounds(tree_formatted=trees_formatted[[", i, "]],dCoeff=X[", i, "],dMat=DiffMat_backwards(X[", length(trees) + 1, "]*SEQ^4+X[", length(trees) + 2, "]*SEQ^2+X[", length(trees) + 3, "]*SEQ),bounds=bounds)", sep = "")
        }
        fun_text <- paste(fun_text, ")}", sep = "") # the end parenthesis
        fun <- eval(parse(text = fun_text))
    }
    return(list(fun = fun, ncoeff = ncoeff, npar = npar, par_names = par_names, par_fixed = list(a = a, b = b, c = c, bounds = bounds), trees = trees, traits = traits, Npts = Npts))
}



testbFPK4 <- lnl_FPK_multiclades_same_V_different_sig2(
    trees = trees[good_indices],
    traits = TRAITS[good_indices],
    bounds = bounds,
    a = NULL, b = NULL, c = NULL, Npts = 50
)

fitbFPK4 <- find.mle_FPK_multiple_clades_same_V_different_sig2(
    model = testbFPK4,
    method = "Nelder-Mead",
    init.optim = NULL
)



fitFPK4$aic
fitbFPK4$aic
fitmFPK4$aic

get.landscape.FPK(fit = fitFPK4)
