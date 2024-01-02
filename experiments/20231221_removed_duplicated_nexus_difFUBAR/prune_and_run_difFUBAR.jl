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


cd("/home/patrick/git/computationalPhylogenetics/experiments/20231221_removed_duplicated_nexus_difFUBAR")


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



for node in getleaflist(tree)
    if node.name in nodes_not_in_fasta
        parent = node.parent

        print(node.name)
        sibling = siblings(node)[1]

        if sibling.children != FelNode[]
            sibling.children[1].parent = parent
            sibling.children[2].parent = parent
            parent.children = sibling.children
        else
            node.parent.name = sibling.name # why does node not have children?
            node.parent.branchlength = sibling.branchlength
            node.parent.children = sibling.children
        end
    end
end

pattern = r"\{([^\}]+)\}"
for node in getleaflist(tree)
    node.name = uppercase(node.name)
    match_result = match(pattern, node.name)
    specie = replace(node.name, pattern => s -> "")
    group = match_result.match

    specie = replace(specie, "-" => s -> "_")
    specie = replace(specie, "__" => s -> "_")
    specie = replace(specie, r"_\z" => s -> "") # remove "_" from last position

    print(specie)
    println(group)
    println(node.name)
    node.name = specie * group
end


tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), MolecularEvolution.universal_code, verbosity=1)

GTRmat # fail 


