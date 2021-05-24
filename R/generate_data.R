# TODO make a
# * genotype matrix with 700,000 variants
# * a phenotype matrix for 1,000 samples
# * a phylogenetic tree

# Phenotype: 
#         Antibiotic_resistance
# sample_1	0
# sample_2	0
# sample_3	1
# sample_4	1

#           SNP_1	SNP_2	SNP_3	SNP_4	SNP_5
# sample_1	0	1	1	0	0
# sample_2	0	0	0	1	1
# sample_3	1	0	0	1	0
# sample_4	1	1	1	0	1

library(ape)
library(phytools)
num_samples <- 1000
num_variants <- 100
set.seed(0)
tree <- phytools::midpoint.root(ape::rcoal(num_samples))
if (!is.rooted(tree)) { stop("not rooted") }

set.seed(0)
bm_phenotype <- abs(floor(phytools::fastBM(tree = tree)))
bm_phenotype[bm_phenotype < 0 | bm_phenotype > 1] <- 0
phenotype <- matrix(data = bm_phenotype, nrow = num_samples, ncol = 1)
row.names(phenotype) <- tree$tip.label
colnames(phenotype) <- "abxR"

genotype <- matrix(data = NA, nrow = num_samples, ncol = num_variants)
row.names(genotype) <- tree$tip.label
colnames(genotype) <- paste0("SNP_", 1:num_variants)
prob_vec <- seq(from = 0.01, to = 0.99, by = 0.01)
set.seed(0)
for (i in 1:ncol(genotype)) {
  probability <- sample(x = prob_vec, size = 1, replace = FALSE, prob = rep(1/99, 99))
  genotype[, i] <- rbinom(n = num_samples, size = 1, prob = probability)
}
write.csv(phenotype, file = "../data/phenotype.csv")
write.csv(genotype, file = "../data/genotype.csv")
write.tree(tree, file = "../data/tree.tree")