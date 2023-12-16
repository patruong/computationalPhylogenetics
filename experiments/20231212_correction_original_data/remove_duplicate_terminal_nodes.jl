using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse
using CSV
using DataFrames
using MolecularEvolution, Plots, Phylo


cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231212_correction_original_data/")


fasta = "ParvoVP_removed_duplicates.fasta" #parsed_args["FASTAfile"]
tree_file = "ParvoVP.nex" #parsed_args["treefile"]
analysis_name = "output" #parsed_args["output_dir"] * "/" * "ParvoVP"
seqnames, seqs = read_fasta(fasta)
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
time_difFUBAR_init = @elapsed tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code

length(getleaflist(tree))



# Generate a list of all sequence name prefixes
seqnames_prefix = []
for seqname in seqnames
    res = split(seqname, "_", limit=2)
    push!(seqnames_prefix, res[1])
end

# Generate a list of nodes that are not in the fasta file but are in the tree file
nodes_not_in_fasta = []
for node in getleaflist(tree)
    if node.children == FelNode[]
        node_name_prefix = split(node.name, "_", limit=2)[1] # to handle various caseing issues
        if !(uppercase(node_name_prefix) in seqnames_prefix)
            push!(nodes_not_in_fasta, node.name)
        end
    end
end

length(seqnames)
length(nodes_not_in_fasta)


length(getleaflist(tree))

for node in getleaflist(tree)
    if node.name in nodes_not_in_fasta
        new_parent = node.parent
        new_name = siblings(node)[1].name
        new_branchlength = siblings(node)[1].branchlength
        new_parent.children = FelNode[]
        new_parent.name = new_name
        new_parent.branchlength = new_branchlength
    end
end



length(getleaflist(tree))




#########################
rest = []
for node in getleaflist(tree)
    if node.name in rest_nodes
        println(node.name)
        push!(rest, node)
    end
end















siblings(rest[1])


rest[1] # find the sibling of this
rest[1].name
n = "MW811187_1_Feline_panleukopenia_virus_strain_FPV-S{G1}"

specific = []
for node in getleaflist(tree)
    if node.name == n
        push!(specific, node)
    end
end

siblings(specific[1])[1].name
siblings(specific[1])[1].name in nodes_not_in_fasta
specific[1].name in nodes_not_in_fasta

println()







pruned_tree = []
function traverse_tree(tree, prune_list)
    for node in tree.children
        if node.name in prune_list
            new_parent = node.parent
            sibling_node = siblings(node)[1]
            #new_parent = sibling_node.parent.parent
            new_parent.children = sibling_node.children
            new_parent.branchlength = sibling_node.branchlength
            new_parent.name = sibling_node.name
            new_parent.nodeindex = sibling_node.nodeindex
            new_parent.seqindex = sibling_node.seqindex
            push!(pruned_tree, new_parent)
        end
        traverse_tree(node)
    end
    return pruned_tree
end
traverse_tree(tree, nodes_not_in_fasta)

rest_nodes = []
for node in getleaflist(tree)
    if node.name in nodes_not_in_fasta
        push!(rest_nodes, node.name)
    end
end

traverse_tree(tree, rest_nodes)


for node in getleaflist(tree)
    if node.name in nodes_not_in_fasta
        println(node.name)
    end
end

unique(node_names)


println("")

println("")

#tree.children[1].children[1]
#tree = tree.children[1]
#tree = siblings(tree.children[1])[1]
#tree.children[1].children[1].children[1].children


# Good

#nodes_not_in_fasta
#function prune_non_existing_terminal_node(tree, nodes_not_in_fasta)
#    node_name = tree.name
#
#    for (child_i, child_node) in enumerate(tree.children)
#        node_name = child_node.name
#        sibling_node = siblings(child_node)[1]
#        sibling_node_name = sibling_node.name
#
#        if node_name in nodes_not_in_fasta
#            if sibling_node_name in nodes_not_in_fasta
#                tree.children = FelNode[] # Remove both leaves
#            else
#                tree = sibling_node
#            end
#        else
#            child_node = prune_non_existing_terminal_node(child_node, nodes_not_in_fasta)
#        end
#    end
#
#    return tree  # Return the modified tree
#end
#
#new_tree = prune_non_existing_terminal_node(tree, nodes_not_in_fasta)
#

### new attempt


getnodelist(tree)

nodes_to_prune = [] # in felsenstein, when they match set the parent to it's sibling.
nodes_to_prune_parents = []
for node in getnodelist(tree)
    if node.name in nodes_not_in_fasta
        push!(nodes_to_prune, node)
        push!(nodes_to_prune_parents, node.parent)
    end
end

length(nodes_to_prune)
length(nodes_to_prune_parents)
# two list one with sibling to prune and one without sibling to prune

remove_both_children_nodes = []
for node in nodes_to_prune
    if siblings(node)[1].name in nodes_not_in_fasta
        push!(remove_both_children_nodes, node.parent)
    end
end

length(remove_both_children_nodes)

function check_children_in_list(tree)
    for child in tree.children
        if child in nodes_not_in_fasta
            return true
        else
            return false
        end
    end
end




#pruned_tree = []
#function traverse_tree(node)
#    for childnode in node.children
#        if childnode in nodes_to_prune_parents
#            println(childnode.name)
#            childnode.children = FelNode[]
#            push!(pruned_tree, childnode)
#        end
#        traverse_tree(childnode)
#    end
#    return pruned_tree
#end
#prune_test = traverse_tree(tree)
#prune_test[33].name
#length(getleaflist(tree))
## this actually worked :/
#


#img = tree_draw(tree)
#savefig_tweakSVG("piechart_tree.svg",img)
#img |> SVG("_tagged_input_tree.svg", 15cm, (3 + length(getleaflist(tree)) / 5)cm)

check_name = "MW811187_1_Feline_panleukopenia_virus_strain_FPV-S{G1}"
#would be good to print the tree
pruned_tree = []
function traverse_tree(node)
    for childnode in node.children
        if childnode.name == check_name
            println("QQQQQQQQQQQQQQQQQQQQQQQQWE ARE HERE")
            println("QQQQQQQQQQQQQQQQQQQEND")
        end

        if childnode in nodes_to_prune_parents
            parent = childnode.parent
            if childnode.children[1].name in nodes_not_in_fasta
                substitution_node = childnode.children[2]
            elseif childnode.children[2].name in nodes_not_in_fasta
                substitution_node = childnode.children[1]
            else
                println("ERROR")
                println(childnode.children[1].name)
                println(childnode.children[1].name in nodes_not_in_fasta)
                println(childnode.children[2].name)
                println(childnode.children[2].name in nodes_not_in_fasta)
                print("CHECK ABOVE")
            end

            childnode.parent = parent
            childnode.children = substitution_node.children
            childnode.branchlength = substitution_node.branchlength
            childnode.name = substitution_node.name
            childnode.nodeindex = substitution_node.nodeindex
            childnode.seqindex = substitution_node.seqindex
            push!(pruned_tree, childnode)
        end
        traverse_tree(childnode)
    end
    return pruned_tree
end
prune_test = traverse_tree(tree)



# some nodes just don't get id'd

length(prune_test)
length(getleaflist(tree)) #we need to set the parent node to the child node

for node in getleaflist(tree)
    if node.name in nodes_not_in_fasta
        println(node.name)
        if node.name == ""
            println(node.children)
        end
    end
end
# why do we still have nodes_not_in_fasta
# we have 12 nodes in fasta still



for node in getleaflist(tree)
    if node.name == "{G2}"
        println(node.name)
        println(node.children)
        println(node)
    end
end

tree




getleaflist(tree)



function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
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

pure_clades = []
function traverse_tree_to_check_for_pure_clades(pure_clades, tree)
    for child in tree.children
        if check_purity_from_node_and_forward_in_time(child)
            #print(child.name)
            #println(" PURE")
            # add this node, as everything under this node is pure
            push!(pure_clades, child)
        else
            #print(child.name)
            #println(" not pure") # keep recursion to se if we go forward in time if we can find pure clades
            traverse_tree_to_check_for_pure_clades(pure_clades, child)
        end
    end
    return pure_clades
end

subtree_tops = traverse_tree_to_check_for_pure_clades(pure_clades, tree)

subtree_tops
subtree_tops

cached_values = Dict()
for x in subtree_tops
    #x = get_node_by_name(x, tree)
    cached_values[x] = Dict()
    parent = x.parent
    x.parent = nothing
    for cp in codon_param_vec
        alpha = cp[1]
        omegas = cp[2:end]

        relevant_omega = omegas[model_ind(x.name, tags)] #here we need to change to child omega

        if haskey(cached_values[x], (alpha, relevant_omega))
            continue
        end

        models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)
        felsenstein!(x, models)
        cached_values[x][(alpha, relevant_omega)] = deepcopy(x.message)
    end
    x.parent = parent
    x.children = FelNode[]
end
# Here i am modifying the specific objects in the subtree_tops list


















function prune_tree(tree)
    for child in tree.children
        if check_children_in_list(child) == true
            # children of this node has stuff to be pruned
            child_1 = child.children[1]
            child_2 = child.children[2]
            if child_1 in nodes_not_in_fasta #remove child_1
                child.children = FelNode[]
                child.parent = child_2
            elseif child_2 in nodes_not_in_fasta #remove child_2
                child.children = FelNode[]
                child.parent = child_1
            else
                println("ERROR")
            end
        else
            prune_tree(child)
        end
    end
end

prune_tree(tree)


for node in getnodelist(tree)
    if node.name in nodes_not_in_fasta
        println(node.name)
        node.parent = FelNode()
    end
end









####
@time tree_, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

tree == tree
tree_ == tree

tree_no_prune = deepcopy(tree)
tree_prune_1 = deepcopy(tree)
tree_prune_2 = deepcopy(tree)
tree_prune_3 = deepcopy(tree)
tree_prune_4 = deepcopy(tree)
tree_prune_5 = deepcopy(tree)

times = []

time_difFUBAR_grid = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 50s
push!(times, time_difFUBAR_grid)

time_difFUBAR_grid_prune_max = @elapsed con_lik_matrix_prune_1, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_1, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s
push!(times, time_difFUBAR_grid_prune_max)

time_difFUBAR_grid_prune_patrick = @elapsed con_lik_matrix_prune_2, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
push!(times, time_difFUBAR_grid_prune_patrick)

time_difFUBAR_grid_prune_patrick_max = @elapsed con_lik_matrix_prune_3, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_3(tree_prune_3, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s
push!(times, time_difFUBAR_grid_prune_patrick_max)

time_difFUBAR_grid_prune_patrick_max_child = @elapsed con_lik_matrix_prune_4, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_4(tree_prune_4, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
push!(times, time_difFUBAR_grid_prune_patrick_max_child)

time_difFUBAR_grid_prune_final = @elapsed con_lik_matrix_prune_5, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_final(tree_prune_5, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
push!(times, time_difFUBAR_grid_prune_final)
