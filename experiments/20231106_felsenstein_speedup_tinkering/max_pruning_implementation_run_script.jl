using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo


# Ace2_tiny
cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
analysis_name = "Ace2"
seqnames, seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

# Simulated data

# HyPhy Sim
cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
file_dir = "/home/ptruong/git/computationalPhylogenetics/contrastFEL_data/omnibus/" # original data
file_dir = "/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering/" # truncated
fasta = file_dir * "sims.1.settings.replicate.1"
tree_file = file_dir * "sims.1.trunc.mod2.nwk"

# Own Sim
# cd("/home/ptruong/git/computationalPhylogenetics/SimulateMG94")
# file_dir = "/home/ptruong/git/computationalPhylogenetics/SimulateMG94/data"
# fasta = file_dir * "sims.1.settings.replicate.1"
# tree_file = file_dir * "sims.1.nwk"

seqnames, seqs = import_hyphy_simulated_FASTA(fasta)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)

comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

# Printout for information
group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

# group_name printout
group_meta_data_file_name = "group_meta_data.txt"
group_meta_data_file = open(group_meta_data_file_name, "w")
for element in group_name
    println(group_meta_data_file, "Groups")
    println(group_meta_data_file, element)
end
println(group_meta_data_file, "")

# pairwise comparison printout
for element in pairwise_comparisons_names
    println(group_meta_data_file, "Comparisons")
    println(group_meta_data_file, element)
end
close(group_meta_data_file)

tag_pos = [1, 2]
tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
analysis_name = "/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering/"
analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)

#difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name;
#    pos_thresh=0.95, iters=2500, verbosity=1,
#    exports=true, code=MolecularEvolution.universal_code)

treestring

@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, analysis_tags, analysis_tag_colors, exports=false, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s




#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s
@time con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s
@time con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=false)


### tree pruning
# if child are pure it is pure
function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end

# Function to check if everything below node is pure.

# getnodelist and get a specific node from nodelist
function get_node_by_name(query_node::String, tree::FelNode)
    for node in getnodelist(tree)
        if node.name == query_node
            return node
        end
    end
end

function check_purity_from_node_and_forward_in_time(tree)
    # collect all nodes under and see if they are the same group
    # if return true that means everything under the tree is pure (the tree node can be differnt though)
    node_groups = []
    for node in getnodelist(tree)
        if tree.name != node.name # we only check nodes under the input tree or subtree
            push!(node_groups, model_ind(node.name, tags))
        end
    end
    if length(unique(node_groups)) == 1
        return true
    else
        return false
    end
end

node_test = get_node_by_name("N363{TEST}", tree)
check_purity_from_node_and_forward_in_time(node) #cool this works

function traverse_tree_to_check_for_pure_clades(tree)
    for child in tree.children
        if check_purity_from_node_and_forward_in_time(child)
            print(child.name)
            println(" PURE")
            # add this node, as everything under this node is pure
        else
            print(child.name)
            println(" not pure") # keep recursion to se if we go forward in time if we can find pure clades
            traverse_node_to_check_for_pure_clades(child)
        end
    end
end

traverse_node_to_check_for_pure_clades(tree)








