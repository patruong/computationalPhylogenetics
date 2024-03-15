
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution
using CodonMolecularEvolution

cd("/home/patrick/git/computationalPhylogenetics/")

#function import_grouped_label_tree(tree_file)
# Takes a Newick tree file and return Newick tree, Newick tree with replaced tags, group tags, original tags, and randomly generated colours for each tag
#tree = read_newick_tree(tree_file)
#treestring = newick(tree)
#treestring_group_labeled, group_tags, original_tags = CodonMolecularEvolution.replace_newick_tags(treestring)
#tag_colors = generate_hex_colors(length(original_tags))
#    return treestring_group_labeled, treestring, group_tags, original_tags, tag_colors
#end


function get_unique_tags(json_file)
    json_data = JSON.parsefile(json_file)
    treestring = json_data["input"]["trees"]["0"]
    tags = json_data["tested"]["0"]
    return unique(values(tags))
end

function read_in_tree_from_hyphy_result_json(json_file)
    json_data = JSON.parsefile(json_file)
    treestring = json_data["input"]["trees"]["0"]
    tags = json_data["tested"]["0"]
    return treestring, tags
end

function tag_hyphy_tree(treestring, tags, group1, group2)
    pattern = r"[\(\),:]([^,:]+)(?=:)"
    matches = eachmatch(pattern, treestring)
    modified_treestring = treestring


    function get_last_match_from_regex(matches)
        last_match = nothing
        for match in matches
            last_match = match.match
        end
        return last_match
    end

    for match in matches
        original_string = match.captures[1]
        regex = r"[A-Za-z0-9_]+"
        pruned_match = eachmatch(regex, original_string)
        key = get_last_match_from_regex(pruned_match)
        #println(key)
        try
            tag = tags[key]
            if tag in [group1, group2]
                replacement = original_string * "{$tag}" # Note this needs to be original_string, otherwise we interpret "Node22" as "Node2"*"2", the original string contains some crap before the string which makes it kinda unique.
            else
                replacement = original_string
            end
            modified_treestring = replace(modified_treestring, original_string => replacement)
        catch
            println("No key $key")
        end
    end

    return modified_treestring
end


# function to read and tag trees without branchlenght
function tag_tree_without_branchlength(treestring, group1, group2, branchlength=0)
    """
    Do not use on a tree with branchlength as this will remove branchlength information
    """
    tree = gettreefromnewick(treestring, FelNode)
    MolecularEvolution.binarize!(tree)
    MolecularEvolution.ladderize!(tree)
    for n in getnodelist(tree)
        n.name = replace(n.name, "'" => "")
        n.branchlength = rand()  # Appends a colon and a random float to the node name for id 
    end
    treestring = newick(tree)
    ## special tagger
    # 2 group + background
    #fasta_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex"
    #json_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex.FEL.json"
    # seqnames, seqs = read_fasta(fasta_file)

    tagged_treestring = tag_hyphy_tree(treestring, tags, group1, group2)
    tree = gettreefromnewick(tagged_treestring, FelNode)
    for n in getnodelist(tree)
        n.name = replace(n.name, "'" => "")
        n.branchlength = branchlength  # set fix value for branchlength
    end
    treestring = newick(tree)
    return treestring
end

function write_string_to_file(filename, str)
    file = open(filename, "w")
    write(file, str)
    close(file)
end


#########################
# 2 group + background  #
# HIV envelope          #
#########################

fasta_file = "data/contrastFEL_empirical_data/hiv-1_envelope/HIV-env.fasta"
json_file = "data/contrastFEL_empirical_data/hiv-1_envelope/HIV-env.fasta.FEL.json"
seqnames, seqs = read_fasta(fasta_file)
treestring, tags = read_in_tree_from_hyphy_result_json(json_file)
get_unique_tags(json_file)
treestring = tag_hyphy_tree(treestring, tags, "HSX", "MSM")
treestring_group_labeled, group_tags, tags = CodonMolecularEvolution.replace_newick_tags(treestring)
tag_colors = CodonMolecularEvolution.generate_hex_colors(length(tags))


# write tree for contrastFEL
tree = gettreefromnewick(treestring, FelNode)
MolecularEvolution.write_nexus("experiments/20240308_contrastFEL_empirical/data/hiv-1_envelope.nex", tree)
write_string_to_file("experiments/20240308_contrastFEL_empirical/data/hiv-1_envelope.nwk", treestring)


tags
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "hiv-1_envelope/analysis"
df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)



# 1 group + background
fasta_file = "data/contrastFEL_empirical_data/epidermal_leaf_trichomes_BRT/BLTFull.fas"
json_file = "data/contrastFEL_empirical_data/epidermal_leaf_trichomes_BRT/BLTFull.fas.FEL.json"
seqnames, seqs = read_fasta(fasta_file)
treestring, tags = read_in_tree_from_hyphy_result_json(json_file)
get_unique_tags(json_file)
treestring = tag_hyphy_tree(treestring, tags, "Physaria", "background")
treestring_group_labeled, group_tags, tags = CodonMolecularEvolution.replace_newick_tags(treestring)
tag_colors = CodonMolecularEvolution.generate_hex_colors(length(tags))

# write tree for contrastFEL
tree = gettreefromnewick(treestring, FelNode)
MolecularEvolution.write_nexus("experiments/20240308_contrastFEL_empirical/data/epidermal_leaf_trichomes_BRT.nex", tree)
write_string_to_file("experiments/20240308_contrastFEL_empirical/data/epidermal_leaf_trichomes_BRT.nwk", treestring)



tags
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "epidermal_leaf/analysis"
df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)

# 2 group + background
fasta_file = "data/contrastFEL_empirical_data/rubisco_C3_vs_C4/rbcl.nex"
json_file = "data/contrastFEL_empirical_data/rubisco_C3_vs_C4/rbcl.nex.FEL.json"
seqnames, seqs = read_fasta(fasta_file)
treestring, tags = read_in_tree_from_hyphy_result_json(json_file)
get_unique_tags(json_file)
treestring = tag_hyphy_tree(treestring, tags, "C3", "C4")
treestring_group_labeled, group_tags, tags = CodonMolecularEvolution.replace_newick_tags(treestring)
tag_colors = CodonMolecularEvolution.generate_hex_colors(length(tags))


# write tree for contrastFEL
tree = gettreefromnewick(treestring, FelNode)
MolecularEvolution.write_nexus("experiments/20240308_contrastFEL_empirical/data/rubisco_C3_vs_C4.nex", tree)
write_string_to_file("experiments/20240308_contrastFEL_empirical/data/rubisco_C3_vs_C4.nwk", treestring)



tags
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "rubisco/analysis"
df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)



## special tagger
# 2 group + background
fasta_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex"
json_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex.FEL.json"
seqnames, seqs = read_fasta(fasta_file)



treestring, tags = read_in_tree_from_hyphy_result_json(json_file)
get_unique_tags(json_file)
treestring = tag_tree_without_branchlength(treestring, "NAIVE", "TREATED", 1)


## EDIT

group1 = "NAIVE"
group2 = "TREATED"

tree = gettreefromnewick(treestring, FelNode)
MolecularEvolution.binarize!(tree)
MolecularEvolution.ladderize!(tree)
for n in getnodelist(tree)
    n.name = replace(n.name, "'" => "")
    n.branchlength = rand()  # Appends a colon and a random float to the node name for id 
end
treestring = newick(tree)

tagged_treestring = tag_hyphy_tree(treestring, tags, group1, group2)
tree = gettreefromnewick(tagged_treestring, FelNode)

tagged_treestring

branchlength = 1
for n in getnodelist(tree)
    n.name = replace(n.name, "'" => "")
    n.branchlength = branchlength  # set fix value for branchlength
end

tree.children[1].branchlength

treestring = newick(tree)





















write_string_to_file("experiments/20240308_contrastFEL_empirical/data/HIV_RT_branchlength_1.nwk", treestring)

# 4 group + background # This is a bit harder to tag
fasta_file = "data/contrastFEL_empirical_data/cytochrome_B_of_Haemosporidians/Cytb.fasta"
json_file = "data/contrastFEL_empirical_data/cytochrome_B_of_Haemosporidians/Cytb.fasta.FEL.json"
seqnames, seqs = read_fasta(fasta_file)

treestring, tags = read_in_tree_from_hyphy_result_json(json_file)

cleaned_seqs = []
for i in 1:length(seqs)
    seq = replace(seqs[1], " " => "")
    push!(cleaned_seqs, seq)
end
seqs = cleaned_seqs
seqs = [string(x) for x in seqs if isa(x, String)]

typeof(seqs)
get_unique_tags(json_file)

treestring
treestring = tag_tree_without_branchlength(treestring, "mammals", "Leucocytozoon", 1)
treestring = tag_tree_without_branchlength(treestring, "mammals", "background", 1)
treestring = tag_tree_without_branchlength(treestring, "mammals", "birds", 1)
treestring = tag_tree_without_branchlength(treestring, "mammals", "Haemoproteidae", 1)
treestring_group_labeled, group_tags, tags = CodonMolecularEvolution.replace_newick_tags(treestring)
tag_colors = CodonMolecularEvolution.generate_hex_colors(length(tags))


tags

# write tree for contrastFEL
tree = gettreefromnewick(treestring, FelNode)
MolecularEvolution.write_nexus("experiments/20240308_contrastFEL_empirical/data/cytb.nex", tree)
write_string_to_file("experiments/20240308_contrastFEL_empirical/data/cytb_branchlength_1.nwk", treestring)


################
# Run difFUBAR #
################



tags
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "cytochrome_B/mammals_leucocytozoon/analysis"
df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)
df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)




