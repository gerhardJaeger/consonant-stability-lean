library(tidyverse)
library(ape)
library(phytools)


d <- read_csv("../data/data_pruned.csv")
tree <- read.nexus("../data/tree_pruned.nex")


ie_tree <- read.nexus("../families/Indo-European/phylogeny/heggartyetal2023.tree")
ie_conversion <- read_csv("../families/Indo-European/phylogeny/IE_conversion.csv") %>%
    distinct(glottocode, .keep_all = TRUE) %>%
    filter(glottocode %in% d$glottocode)


an_tree <- read.nexus("../families/Austronesian/phylogeny/grayetal2009.summary.trees")
an_conversion <- read_csv("../families/Austronesian/phylogeny/austronesian_conversion.csv") 


bantu_tree <- read.nexus("../families/Bantu/phylogeny/ctmc4g-asc.ucln.yule.ba-sp.tree")
bantu_tree$edge.length <- bantu_tree$edge.length / 1000
bantu_conversion <- read_csv("../families/Bantu/phylogeny/Bantu_Glottocodes_Koile.csv")


pn_tree <- read.nexus("../families//Pama-Nyungan/phylogeny/bouckaertetal2018.txt")
pn_tree$edge.length <- pn_tree$edge.length / 1000
pn_conversion <- read_csv("../families/Pama-Nyungan/phylogeny/PN_conversion.csv")

st_tree <- read.nexus("../families/Sino-Tibetan/phylogeny/Zhangetal2019.MCC.tree")
st_tree$edge.length <- st_tree$edge.length / 1000
st_conversion <- read_csv("../families/Sino-Tibetan/phylogeny/Zhang_ST_conversion.csv")


tg_tree <- read.nexus("../families/Tupi-Guarani/phylogeny/gerardi2023.full.relaxed.mcc.nex")
tg_tree$edge.length <- tg_tree$edge.length / 1000
tg_conversion <- read_csv("../families/Tupi-Guarani/phylogeny/TG_conversion.csv")

turkic_tree <- read.nexus("../families/Turkic/phylogeny/savelyevetal2020_turkic_AnnotatedTree.trees")
turkic_conversion <- read_csv("../families/Turkic/phylogeny/Turkic_conversion.csv")


uralic_tree <- read.nexus("../families/Uralic/phylogeny/uralic-no_constraint.mcc.nex")
uralic_tree$edge.length <- uralic_tree$edge.length / 1000
uralic_conversion <- read_csv("../families/Uralic/phylogeny/uralic_conversion.csv")

##

ua_tree <- read.nexus("../families/Uto-Aztecan/phylogeny/greenhill2023-covarion-relaxed.mcct.trees")
ua_tree$edge.length <- ua_tree$edge.length / 1000
ua_conversion <- read_csv("../families/Uto-Aztecan/phylogeny/UA_conversion.csv") %>%
    filter(glottocode != "kiow1266" )


##
replace_subtree <- function(tree, fam_tree, fam_conversion) {
    fam_conversion <- fam_conversion %>% 
        distinct(glottocode, .keep_all = TRUE) %>%
        filter(glottocode %in% tree$tip.label)
    fam_tree <- drop.tip(fam_tree, setdiff(fam_tree$tip.label, fam_conversion$name))
    fam_tree$tip.label <- fam_conversion[match(
        fam_tree$tip.label,
        fam_conversion$name
    ), "glottocode"] %>% pull
    fam_root <- getMRCA(tree, fam_conversion$glottocode)
    pre_fam <- getParent(tree, fam_root)
    pre_fam_depth <- node.depth.edgelength(tree)[pre_fam]
    age_proto_world <- max(node.depth.edgelength(tree))
    fam_age <- node.depth.edgelength(fam_tree) %>% max
    pruned_tree <- bind.tip(
        tree,
        "Proto-Family",
        edge.length = age_proto_world - pre_fam_depth - fam_age,
        where = pre_fam
    )
    fam_root <- getMRCA(pruned_tree, fam_conversion$glottocode)
    fam_leaves <- intersect(
        getDescendants(pruned_tree, fam_root),
        1:length(pruned_tree$tip.label)
    )
    pruned_tree <- drop.tip(pruned_tree, fam_leaves)
    new_tree <- bind.tree(
        pruned_tree,
        fam_tree,
        where = which(pruned_tree$tip.label == "Proto-Family")
    )
    return(new_tree)
}

##

tree <- tree %>% 
    replace_subtree(ie_tree, ie_conversion) %>%
    replace_subtree(an_tree, an_conversion) %>%
    replace_subtree(bantu_tree, bantu_conversion) %>%
    replace_subtree(pn_tree, pn_conversion) %>%
    replace_subtree(st_tree, st_conversion) %>%
    replace_subtree(tg_tree, tg_conversion) %>%
    replace_subtree(turkic_tree, turkic_conversion) %>%
    replace_subtree(uralic_tree, uralic_conversion) %>%
    replace_subtree(ua_tree, ua_conversion)
    
tree <- collapse.singles(tree)
tree <- multi2di(tree)

write.nexus(tree, file = "tree_families_replaced.nex")

