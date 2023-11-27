using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo




using Pkg
Pkg.activate("environments/CodonMolecularEvolution")
#using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
#using MolecularEvolution, FASTX, CodonMolecularEvolution

#using Pkg
#Pkg.add("Revise")
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution

# HyPhy Sim
cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
file_dir = "/home/ptruong/git/computationalPhylogenetics/contrastFEL_data/omnibus/" # original data
file_dir = "/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering/" # truncated
fasta = file_dir * "sims.1.settings.replicate.1"
tree_file = file_dir * "sims.1.trunc.mod.nwk"

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
analysis_name = "/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering"
analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)

#difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name;
#    pos_thresh=0.95, iters=2500, verbosity=1,
#    exports=true, code=MolecularEvolution.universal_code)

treestring

@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, analysis_tags, analysis_tag_colors, exports=false, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s


#######################
# Work with tree here #
#######################

img = tree_draw(tree)
img |> PDF("imgout_sim.pdf", 10cm, 10cm)



treestring

tree.children[1].name

# figure out how to traverse down a subtree.
# figure out how to get a pure clade.
x = tree.children[2].children[1].children[1].children[1].children[1].children[1].children[1].children[1]
x_empty = tree.children[2].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children


function extract_group(str)
    curly_bracket_pattern = r"{([^}]+)}"
    extracted_group = match(curly_bracket_pattern, str)
    if isnothing(extracted_group)
        return "" # return "background" 
    else
        return extracted_group[1]
    end
end

# recursion step 1 of 8
tree.children[2]s

# left child
l_child = tree.children[2].children[1].children[1].children[1].children[1].children[2].children[1].children[1]
# right child
r_child = tree.children[2].children[1].children[1].children[1].children[1].children[2].children[1].children[2]


r_child.parent
l_group = extract_group(l_child.name)
r_group = extract_group(r_child.name)

l_group == r_group

r_child.pare


fieldnames(typeof(x))
isempty(fieldnames(typeof(x)))
isempty(fieldnames(typeof(x_empty)))
isempty(x_empty)

tree

tree.children[2].children[1].children[1].children[1].children[1].children[2].children[1].children[1].name
tree.children[2].children[1].children[1].children[1].children[1].children[2].children[1].children[2].name



tree.children[2].children[1].children[1].children[1].children[1].children[1].children[1].children[2].name







img = tree_draw(tree)
img |> PDF("imgout_sim_2group.pdf", 10cm, 70cm)


#####################
# rest of the code ##
#####################



con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=false)



##################### #########
# Working with the treestring #
###############################


treestring
curly_bracket_pattern = r"{([^}]+)}"
extracted_group = match(curly_bracket_pattern, treestring)

groups = findall(curly_bracket_pattern, treestring)

treestring[groups[1]]
treestring

(root
    (level1
        (level2
            (leve3
            )
        )
    )
)

count the brackets and assign levels, and then do something to make it identify subtrees

or can we begin from leaves? if left and right is pure clades and we merge against something that is pure it is still pure...
