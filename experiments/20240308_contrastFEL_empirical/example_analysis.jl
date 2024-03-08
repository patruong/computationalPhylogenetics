
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

function read_in_tree_from_hyphy_result_json(json_file, group1, group2, tags_dict)
    json_data = JSON.parsefile(json_file)
    treestring = json_data["input"]["trees"]["0"]

    for (key, value) in tags_dict
        if value in [group1, group2]
            # Improved pattern to match key before any separator or colon
            pattern = r"(?<=[\(\),])($key):"
            matches = eachmatch(pattern, treestring)
            for match in matches
                start_index = match.offset
                end_index = start_index + length(match.match) - 2  # Adjust for colon
                tagged_string = treestring[start_index:end_index] * "{" * value * "}"
                treestring = replace(replace(treestring, match.match => ""),
                                   treestring[start_index:end_index] => tagged_string)
            end
        end
    end
    return treestring
end




json_data = JSON.parsefile(json_file)
treestring = json_data["input"]["trees"]["0"]




treestring
tags = json_data["tested"]["0"]



pattern = r"[\(\),:]([^,:]+)(?=:)"



matches = eachmatch(pattern, treestring)


modified_treestring = treestring


# Loop through the matches and convert them using the dictionary
for match in matches
    original_string = match.captures[1]
    extra_prune_pattern = r"\(?([^(),:]+)"

    pruned_match= eachmatch(extra_prune_pattern, original_string)
    key = first(pruned_match).captures[1]


    println(original_string)
    #println(key)
    #replacement_value = get(tags, original_string, "test")
    #print(replacement_value)
    #break
    #replacement = original_string * "{$replacement_value}"
    #modified_treestring = replace(modified_treestring, original_string => replacement)
end


test_string = "((((((((HSX_2_CON_9029_1998"
test_string = "((HSX_2_CON_9029_1998"
test_string = ")HSX_2_CON_9029_1998"
test_string = "0.0004927)Node228" # this does not work yet... make this work so i can parse tag the tree to run difFUBAR and contrastFEL
matches = eachmatch(test_pattern, test_string)
original_string = first(matches).captures[1]
println(original_string)

get(tags, ", "test")
.captures[1]

modified_treestring

tags
tags


# 2 group + background
fasta_file = "data/contrastFEL_empirical_data/hiv-1_envelope/HIV-env.fasta"
json_file = "data/contrastFEL_empirical_data/hiv-1_envelope/HIV-env.fasta.FEL.json"
get_unique_tags(json_file)
read_in_tree_from_hyphy_result_json(json_file, "HSX", "MSM")


# 4 group + background
fasta_file = "data/contrastFEL_empirical_data/cytochrome_B_of_Haemosporidians/Cytb.fasta"
json_file = "data/contrastFEL_empirical_data/cytochrome_B_of_Haemosporidians/Cytb.fasta.FEL.json"
get_unique_tags(json_file)

# 1 group + background
fasta_file = "data/contrastFEL_empirical_data/epidermal_leaf_trichomes_BRT/BLTFull.fas"
json_file = "data/contrastFEL_empirical_data/epidermal_leaf_trichomes_BRT/BLTFull.fas.FEL.json"
get_unique_tags(json_file)

# 2 group + background
nex_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex"
json_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex.FEL.json"
get_unique_tags(json_file)

# 2 group + background
nex_file = "data/contrastFEL_empirical_data/rubisco_C3_vs_C4/rbcl.nex"
json_file = "data/contrastFEL_empirical_data/rubisco_C3_vs_C4/rbcl.nex.FEL.json"
get_unique_tags(json_file)







group1 = "HSX"
group2 = "background"

json_data = JSON.parsefile(json_file)
treestring = json_data["input"]["trees"]["0"]
tags = json_data["tested"]["0"]
unique(values(tags))



treestring

treestring = read_in_tree_from_hyphy_result_json(json_file)



tags
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95

seqnames, seqs = read_fasta(fasta_file)
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
treestring_group_labeled, treestring, group_tags, tags, tag_colors = import_grouped_label_tree(tree_file)



df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)


