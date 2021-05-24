library(ape)

phenotype <- read.csv(file = "../data/phenotype.csv", row.names = 1)
tree <- read.tree(file = "../data/tree.tree")

plot(tree, show.tip.label = FALSE)
tiplabels(pch = 15, tip = (1:Ntip(tree))[phenotype == 1], col = "red", cex = 0.1)
tiplabels(pch = 15, tip = (1:Ntip(tree))[phenotype == 0], col = "blue", cex = 0.1)
