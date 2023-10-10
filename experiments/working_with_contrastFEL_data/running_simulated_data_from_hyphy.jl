using Pkg

Pkg.activate("environments/CodonMolecularEvolution")
Pkg.status()
cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution

#cd("/home/ptruong/git/computationalPhylogenetics/contrastFEL_data/") 
cd("/home/ptruong/git/computationalPhylogenetics/")


function import_hyphy_simulated_FASTA(file_name)
    # Import data generated using https://github.com/veg/hyphy-analyses/tree/master/SimulateMG94
    sequence_names = []
    sequences = []

    open(file_name, "r") do file
        current_sequence = ""

        for line in eachline(file)
            if length(line) == 0
                break
            else
                if line[1] == '>'
                    # Store the current sequence and reset for the next entry
                    if !isempty(current_sequence)
                        push!(sequences, current_sequence)
                    end
                    push!(sequence_names, line[2:end])
                    current_sequence = ""
                else
                    # Append the line to the current sequence
                    current_sequence *= line
                end
            end
        end

        # Store the last sequence after the loop ends (if any)
        if !isempty(current_sequence)
            push!(sequences, current_sequence)
        end
    end
    sequence_names = [string(seq_name) for seq_name in sequence_names]
    sequences = [string(seq) for seq in sequences]

    return sequence_names, sequences
end

function replace_newick_tags(treestr)
    # This function replaces whatever is in {} in a newick tree to G1...Gn

    pattern = r"\{([^}]+)\}"
    unique_tags = Set{String}()
    for match in eachmatch(pattern, treestr)
        println(match)
        push!(unique_tags, match.match)
    end
    unique_tags = collect(unique_tags)
    num_unique_tags = length(unique_tags)
    group_tags = ["{G$index}" for (index, _) in enumerate(unique_tags)]

    if length(unique_tags) != length(group_tags)
        throw(ArgumentError("The number of unique tags must be equal to the number of group tags."))
    end

    tag_mapping = Dict{String,String}()

    for (old_tag, new_tag) in zip(unique_tags, group_tags)
        tag_mapping[old_tag] = "$new_tag"
    end

    for (old_tag, new_tag) in tag_mapping
        treestr = replace(treestr, old_tag => new_tag)
    end

    return treestr, group_tags, unique_tags
end

# Define a function to generate distinct hexadecimal color codes
function generate_hex_colors(num_colors)
    # Generates an arbitrary number of hex-color to match group_tags
    colors = []
    for i in 1:num_colors
        r = rand(0:255)
        g = rand(0:255)
        b = rand(0:255)
        hex_color = "#" * string(r, base=16, pad=2) * string(g, base=16, pad=2) * string(b, base=16, pad=2)
        push!(colors, hex_color)
    end
    colors_string = [string(color) for color in colors]
    return colors_string
end

function import_grouped_label_tree(tree_file)
    # Takes a Newick tree file and return Newick tree, Newick tree with replaced tags, group tags, original tags, and randomly generated colours for each tag
    tree = read_newick_tree(tree_file)
    treestring = newick(tree)
    treestring_group_labeled, group_tags, original_tags = replace_newick_tags(treestring)
    tag_colors = generate_hex_colors(length(original_tags))
    return treestring_group_labeled, treestring, group_tags, original_tags, tag_colors
end

exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "output/test_sim_data"


cd("/home/ptruong/git/computationalPhylogenetics/")
file_name = "SimulateMG94/data/example3/example3.replicate.1.nex"
seqnames, seqs = import_hyphy_simulated_FASTA(file_name)

tree_file = "SimulateMG94/CD2.nwk"
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)

#@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, original_tags, tag_colors, exports=exports, verbosity=verbosity)
@time code = MolecularEvolution.universal_code
@time code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code) #19s
@time con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 73 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)



