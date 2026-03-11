library(tidyverse)
library(ape)
library(treeio)
library(glue)

dir.create("data/asr", recursive = TRUE, showWarnings = FALSE)

# Load main data and tree
d     <- read_csv("data/data_pruned.csv", show_col_types = FALSE)
tree  <- read.nexus("data/tree_families_replaced.nex")

# Load ASR posterior samples (rows = nodes, columns = posterior samples V1, V2, ...)
asr_plosive   <- read_csv("data/asr/asr_plosive.csv",   show_col_types = FALSE)
asr_fricative <- read_csv("data/asr/asr_fricative.csv", show_col_types = FALSE)
asr_affricate <- read_csv("data/asr/asr_affricate.csv", show_col_types = FALSE)
asr_series    <- read_csv("data/asr/asr_series.csv",    show_col_types = FALSE)

# Compute per-node summary statistics across posterior samples
summarise_asr <- function(asr_df) {
  mat <- as.matrix(asr_df)
  tibble(
    node  = seq_len(nrow(mat)),
    mean  = rowMeans(mat),
    lower = apply(mat, 1, quantile, probs = 0.025),
    upper = apply(mat, 1, quantile, probs = 0.975)
  )
}

summary_plosive   <- summarise_asr(asr_plosive)
summary_fricative <- summarise_asr(asr_fricative)
summary_affricate <- summarise_asr(asr_affricate)
summary_series    <- summarise_asr(asr_series)

# Annotate a treedata object and write BEAST nexus
write_annotated_tree <- function(tree, summary_df, out_path) {
  td <- treeio::as.treedata(tree)
  td@data <- summary_df
  treeio::write.beast(td, file = out_path)
}

write_annotated_tree(tree, summary_plosive,   "data/asr/world_plosive.nexus")
write_annotated_tree(tree, summary_fricative, "data/asr/world_fricative.nexus")
write_annotated_tree(tree, summary_affricate, "data/asr/world_affricate.nexus")
write_annotated_tree(tree, summary_series,    "data/asr/world_series.nexus")

# Family definitions and conversion file paths
families <- list(
  `Indo-European`  = "data/families/Indo-European/phylogeny/IE_conversion.csv",
  `Austronesian`   = "data/families/Austronesian/phylogeny/austronesian_conversion.csv",
  `Bantu`          = "data/families/Bantu/phylogeny/Bantu_Glottocodes_Koile.csv",
  `Pama-Nyungan`   = "data/families/Pama-Nyungan/phylogeny/PN_conversion.csv",
  `Sino-Tibetan`   = "data/families/Sino-Tibetan/phylogeny/Zhang_ST_conversion.csv",
  `Tupi-Guarani`   = "data/families/Tupi-Guarani/phylogeny/TG_conversion.csv",
  `Turkic`         = "data/families/Turkic/phylogeny/Turkic_conversion.csv",
  `Uralic`         = "data/families/Uralic/phylogeny/uralic_conversion.csv"
)

trait_asrs <- list(
  plosive   = list(samples = asr_plosive,   summary = summary_plosive),
  fricative = list(samples = asr_fricative, summary = summary_fricative),
  affricate = list(samples = asr_affricate, summary = summary_affricate),
  series    = list(samples = asr_series,    summary = summary_series)
)

n_tips <- length(tree$tip.label)

for (family_name in names(families)) {
  conv_path <- families[[family_name]]
  conv      <- read_csv(conv_path, show_col_types = FALSE)

  # All conversion CSVs have a glottocode column; tree tip labels are glottocodes
  family_glottocodes <- conv$glottocode

  # Keep only glottocodes present in the tree
  tips_to_keep <- intersect(family_glottocodes, tree$tip.label)

  if (length(tips_to_keep) < 2) {
    warning(glue("Family {family_name}: fewer than 2 matching tips found; skipping."))
    next
  }

  # Prune global tree to family tips
  tips_to_drop  <- setdiff(tree$tip.label, tips_to_keep)
  family_tree   <- drop.tip(tree, tips_to_drop)

  # Map pruned tree node indices back to original tree node indices
  # Tips first, then internal nodes
  original_tip_indices   <- match(family_tree$tip.label, tree$tip.label)
  # Internal nodes: identify via mrca on original tree
  family_n_tips    <- length(family_tree$tip.label)
  family_n_nodes   <- family_tree$Nnode
  family_n_total   <- family_n_tips + family_n_nodes

  # Build a mapping from family-tree node index to original-tree node index
  # for tips this is straightforward; for internal nodes use getMRCA
  node_map <- integer(family_n_total)
  node_map[seq_len(family_n_tips)] <- original_tip_indices

  for (fam_node in seq_len(family_n_nodes)) {
    orig_fam_node_idx <- family_n_tips + fam_node
    # descendants of this internal node in the family tree
    desc_tips_fam <- family_tree$tip.label[
      which(family_tree$tip.label %in%
        family_tree$tip.label[unlist(phangorn::Descendants(family_tree, orig_fam_node_idx, type = "tips"))])
    ]
    if (length(desc_tips_fam) >= 2) {
      orig_desc_indices <- match(desc_tips_fam, tree$tip.label)
      node_map[orig_fam_node_idx] <- ape::getMRCA(tree, orig_desc_indices)
    } else if (length(desc_tips_fam) == 1) {
      node_map[orig_fam_node_idx] <- match(desc_tips_fam, tree$tip.label)
    }
  }

  for (trait_name in names(trait_asrs)) {
    asr_samples <- trait_asrs[[trait_name]]$samples

    # Extract rows corresponding to the family nodes (using node_map)
    valid_idx <- node_map[node_map > 0]
    fam_asr   <- asr_samples[valid_idx, , drop = FALSE]

    fam_summary <- tibble(
      node  = seq_len(nrow(fam_asr)),
      mean  = rowMeans(as.matrix(fam_asr)),
      lower = apply(as.matrix(fam_asr), 1, quantile, probs = 0.025),
      upper = apply(as.matrix(fam_asr), 1, quantile, probs = 0.975)
    )

    out_path <- glue("data/asr/{family_name}_{trait_name}.nexus")
    write_annotated_tree(family_tree, fam_summary, out_path)
  }

  message(glue("Done: {family_name}"))
}
