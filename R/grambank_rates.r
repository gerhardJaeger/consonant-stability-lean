library(geiger)
library(tidyverse)
library(pbapply)


##

tree <- read.nexus("../data/grambank/grambank_tree.tre")

d <- read_csv("../data/grambank/values.csv") %>%
    select(Language_ID, Parameter_ID, Value) %>%
    filter(Value != "?") %>%
    filter(Language_ID %in% tree$tip.label) %>%
    mutate(Value = as.numeric(Value)) %>%
    drop_na()

features <- d$Parameter_ID %>% unique

##

min1 <- d %>%
    group_by(Parameter_ID) %>%
    summarize(minvalue = min(Value)) %>%
    filter(minvalue > 0) %>%
    pull(Parameter_ID)


for (f in min1) {
    d[d$Parameter_ID == f, "Value"] <- d[d$Parameter_ID == f, "Value"] - 1
}


d <- d %>%
    filter(Value < 2)



##


characteristic_time <- function(f) {
    d_f <- filter(d, Parameter_ID == f)
    x <- as.vector(d_f$Value) + 1
    names(x) <- as.vector(d_f$Language_ID)
    tree_x <- drop.tip(tree, setdiff(tree$tip.label, names(x)))
    fit_x <- fitDiscrete(tree_x, x, model = "ARD")
    # get the names from fit_x$opt starting with q
    q <- fit_x$opt[grepl("^q", names(fit_x$opt))]
    lambda <- q$q12 + q$q21
    return(1/lambda)
}


##

tc <- pbsapply(features, function(i) characteristic_time(i))

##
parameters <- read_csv("../data/grambank/parameters.csv")

##

tc_estimates <- data.frame(feature = features, tc = tc) %>%
    arrange(desc(tc)) %>%
    left_join(select(parameters, ID, Name), by = c("feature" = "ID"))

write_csv(tc_estimates, "grambank_tc_estimates.csv")



##

d_series <- read_csv("../data/data_pruned.csv")
d_series <- d_series %>%
    filter(glottocode %in% tree$tip.label)
d_series <- d_series[match(tree$tip.label, d_series$glottocode),]


##

series_vars <- c("plosive_markedness_fullness", "fricative_markedness_fullness", "affricate_markedness_fullness", "series_markedness_fullness")


binary_vars <- str_c(series_vars, "_binary")

median_values <- d_series %>%
    select(series_vars) %>%
    summarize(across(everything(), median, na.rm = TRUE))




binarize_series <- function(f) {
    return(ifelse(d_series[[f]] > median_values[[f]], 1, 0))
}


for (i in 1:length(series_vars)) {
    f <- series_vars[i]
    d_series <- d_series %>%
        mutate(
            !!binary_vars[i] := binarize_series(f)
        )
}


##

characteristic_time_series <- function(f) {
    x <- as.vector(d_series[[f]]) + 1
    names(x) <- as.vector(d_series$glottocode)
    tree_x <- drop.tip(tree, setdiff(tree$tip.label, names(x)))
    fit_x <- fitDiscrete(tree_x, x, model = "ARD")
    # get the names from fit_x$opt starting with q
    q <- fit_x$opt[grepl("^q", names(fit_x$opt))]
    lambda <- q$q12 + q$q21
    return(1/lambda)
}


tc_series <- pbsapply(binary_vars, characteristic_time_series)

##


dens <- density(log(2) * tc_estimates$tc)
y_max <- max(dens$y)

ggplot(tc_estimates, aes(x = log(2) * tc)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "grey80", color = "white") +
    geom_density(fill = "blue", alpha = 0.3) +
    labs(
        title = "Estimated Characteristic Time for Grambank Features",
        x = "Estimated Characteristic Time",
        y = "Density"
    ) +
    geom_vline(xintercept = tc_series[1], linetype = "dotdash") +
    geom_vline(xintercept = tc_series[2], linetype = "dashed") +
    geom_vline(xintercept = tc_series[3], linetype = "dotted") +
    geom_vline(xintercept = tc_series[4]) +
    annotate("text", x = tc_series[1], y = y_max * 1.05, label = "plosives", angle = 90, vjust = -0.5, color = "black") +
    annotate("text", x = tc_series[2], y = y_max * 1.05, label = "fricatives", angle = 90, vjust = -0.5, color = "black") +
    annotate("text", x = tc_series[3], y = y_max * 1.05, label = "affricates", angle = 90, vjust = -0.5, color = "black") +
    annotate("text", x = tc_series[4], y = y_max * 1.05, label = "series", angle = 90, vjust = -0.5, color = "black")
