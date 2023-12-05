library("ggtree")
library("treeio")

nwk_tree <- read.tree("tree_newick.nwk")
nwk_tree

nexus_tree <- read.nexus("sim_tree.nex")
nexus_tree

# plotting
#ggtree_object <- ggtree(nexus_tree)
#ggtree_object
#plot(ggtree_object)
