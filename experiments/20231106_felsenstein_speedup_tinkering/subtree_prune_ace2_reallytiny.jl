using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo



cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
analysis_name = "Ace2"
seqnames, seqs = read_fasta("data/Ace2_reallytiny.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/ACE2_reallytiny_tagged.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s
con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=true)

img = tree_draw(tree)
img |> PDF("imgout_ace2_reallytiny.pdf", 10cm, 10cm)


getnodelist(tree)[3].name

function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end

strip_tags_from_name = CodonMolecularEvolution.generate_tag_stripper(tags)

color_dict = Dict(zip(getnodelist(tree), [tag_colors[model_ind(n.name, tags)] for n in getnodelist(tree)]))
label_dict = Dict(zip(getnodelist(tree), [strip_tags_from_name(n.name) for n in getnodelist(tree)]))


# Full tree representation
# Now figure out how to id subtree.
x = tree.name
model_ind(x, tags)

x = tree.children[1].name
model_ind(x, tags)
x = tree.children[1].children[1].name
model_ind(x, tags)
x = tree.children[1].children[2].name
model_ind(x, tags)
x = tree.children[1].children[2].children[1].name
model_ind(x, tags)
x = tree.children[1].children[2].children[2].name
model_ind(x, tags)
x = tree.children[1].children[2].children[2].children[1].name
model_ind(x, tags)
x = tree.children[1].children[2].children[2].children[2].name
model_ind(x, tags)

x = tree.children[2].name
model_ind(x, tags)
x = tree.children[2].children[1].name
model_ind(x, tags)
x = tree.children[2].children[2].name
model_ind(x, tags)
x = tree.children[2].children[2].children[1].name
model_ind(x, tags)
x = tree.children[2].children[2].children[2].name
model_ind(x, tags)
x = tree.children[2].children[2].children[2].children[1].name
model_ind(x, tags)
x = tree.children[2].children[2].children[2].children[2].name
model_ind(x, tags)




function extract_group(str)
    curly_bracket_pattern = r"{([^}]+)}"
    extracted_group = match(curly_bracket_pattern, str)
    if isnothing(extracted_group)
        return "" # return "background" 
    else
        return extracted_group[1]
    end
end

pure_leaf_tag = true

test_str = "testunit{G2}"
test_str2 = "test_nogroup"
extract_group(test_str)
extract_group(test_str2)


x = tree.children[2].children[2].children[2]

current_node_group = nothing # group of the current node
# if child below current node are the same, it is pure. 
# end at check if there is child
left_child_group = extract_group(x.children[1].name)
right_child_group = extract_group(x.children[2].name)
if left_child_group == right_child_group
    println("pure")
end

# Write the whole recursion

current_group = extract_group(tree.name)
extract_group(tree.children[1].name)
extract_group(tree.children[2].name)



### my paper algo

x = tree

current_group = extract_group(x.name)
l_child = x.children[1]
r_child = x.children[2]
l_group = extract_group(l_child.name)
r_group = extract_group(r_child.name)
pure_children = (l_group == r_group)

l_child


# find subtree thing

stack = [(tree, 1, true)]

while length(stack) > 0
    node, ind, first = pop!(stack)
    if !isleafnode(node)
        if first
            push!(stack, (node, ind, false))
            for i = 1:length(node.children)
                push!(stack, (node.children[i], i, true))
            end
        end
        #if !first
        #    for part in partition_list
        #        #Combine child messages into node message.
        #        combine!(
        #            node.message[part],
        #            [mess[part] for mess in node.child_messages],
        #            true,
        #        )
        #    end
        #    if !isroot(node)
        #        #Get the model list for the current branch.
        #        model_list = models(node)
        #        for part in partition_list
        #            backward!(
        #                node.parent.child_messages[ind][part],
        #                node.message[part],
        #                model_list[part],
        #                node,
        #            )
        #        end
        #    end
        #end
    end
end

# understand how to leaf down to id the sub-trees and then it should be pretty straight-forward.
test = [1, 2, 3]
push!(test, 4)