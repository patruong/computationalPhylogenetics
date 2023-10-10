using Pkg

Pkg.activate("environments/CodonMolecularEvolution")
Pkg.status()
cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution

cd("/home/ptruong/git/computationalPhylogenetics/workflow_pipeline/") # how do we fix this to snakemake?

exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "output/results/difFUBAR"

#seqnames, seqs = read_fasta("output/data/flu_sim_MG94_F61.fasta");
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("output/data/flu_sim_after_sample_down.nex")

#seqnames, seqs = read_fasta("output/data/flu_sim_MG94_F61.fasta");
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("output/data/example3.replicate.1.nex")


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

function replace_newick_tags(newick_tree, unique_tags, group_tags)
    # This function replaces whatever is in {} in a newick tree to G1...Gn
    if length(unique_tags) != length(group_tags)
        throw(ArgumentError("The number of unique tags must be equal to the number of group tags."))
    end

    tag_mapping = Dict{String,String}()

    for (old_tag, new_tag) in zip(unique_tags, group_tags)
        tag_mapping[old_tag] = "{$new_tag}"
    end

    for (old_tag, new_tag) in tag_mapping
        newick_tree = replace(newick_tree, old_tag => new_tag)
    end

    return newick_tree
end

function import_hyphy_simulated_data(file_name)
    # Import data generated using https://github.com/veg/hyphy-analyses/tree/master/SimulateMG94
    sequence_names = []
    sequences = []
    treestring = []
    tags = ""
    tag_colors = Dict{String,String}()

    open(file_name, "r") do file
        current_sequence = ""

        for line in eachline(file)
            if length(line) == 0
                nothing
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

    start_substr = "("
    lines = readlines(file_name)
    treeline = lines[findfirst([occursin(start_substr, l) for l in lines])]
    #st = findfirst(start_substr, treeline)[end]
    #treestr = treeline[st+1:end]
    treestr = treeline

    # Regular expression pattern to match tags within curly braces {}
    pattern = r"\{([^}]+)\}"
    unique_tags = Set{String}()
    for match in eachmatch(pattern, treestr)
        println(match)
        push!(unique_tags, match.match)
    end
    unique_tags = collect(unique_tags)
    num_unique_tags = length(unique_tags)
    group_tags = ["G$index" for (index, _) in enumerate(unique_tags)]

    # Generate colors for each group
    num_groups = length(group_tags)
    tag_colors = generate_hex_colors(num_groups)

    treestr = replace_tags(treestr, unique_tags, group_tags) # replace whatever is in {} to G1...Gn

    return sequence_names, sequences, treestr, unique_tags, group_tags, tag_colors

end

file_name = "output/data/example3.replicate.1.nex"
seqnames, seqs, treestring, unique_tags, tags, tag_colors = import_hyphy_simulated_data(file_name)


###
seqnames, seqs = read_fasta("output/data/test_data/Ace2_tiny_tagged.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("output/data/test_data/Ace2_no_background.nex")

#@time df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name);
treestring = "(((XM_027533928_Bos_indicus_x_Bos_taurus{G1},(XM_042974087_Panthera_tigris{G1},(XM_008696415_Ursus_maritimus{G1},(XM_025857612_Callorhinus_ursinus{G1},XM_034996559_Mirounga_leonina{G1}){G1}){G1}){G1}){G1},(((bat_Trachops_cirrhosus{G1},bat_Pteronotus_parnellii{G1}){G1},(bat_Molossus_molossus{G1},bat_Pipistrellus_kuhlii{G1}){G1}){G1},((bat_Rousettus_aegyptiacus{G1},bat_Cynopterus_brachyotis{G1}){G1},(bat_Rhinolophus_sinicus{G1},bat_Hipposideros_pratti{G1}){G1}){G1}){G1}){G1},(((XM_037822651_Choloepus_didactylus{G2},(XM_004866100_Heterocephalus_glaber{G2},(XM_012638731_Propithecus_coquereli{G2},((XM_008988993_Callithrix_jacchus{G2},XM_032285963_Sapajus_apella{G2}){G2},((XM_003261084_Nomascus_leucogenys{G2},((NM_001131132_Pongo_abelii{G2},XM_024240245_Pongo_abelii{G2}){G2}){G2}){G2},((XM_019019204_Gorilla_gorilla_gorilla{G2},(NM_001371415_Homo_sapiens{G2},(XM_008974180_Pan_paniscus{G2},XM_016942979_Pan_troglodytes{G2}){G2}){G2}){G2}){G2}){G2}){G2},((XM_010366065_Rhinopithecus_roxellana{G2},XM_023199053_Piliocolobus_tephrosceles{G2}){G2},((XM_011995533_Mandrillus_leucophaeus{G2},XM_012035808_Cercocebus_atys{G2}){G2},(XM_011735203_Macaca_nemestrina{G2},(NM_001135696_Macaca_mulatta{G2},(XM_005593037_Macaca_fascicularis{G2},XM_015126958_Macaca_mulatta{G2}){G2}){G2}){G2}){G2}){G2}){G2}){G1};"

# There is some data missing from the hyphy simulatiomns

@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
@time code = MolecularEvolution.universal_code
@time code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code) #19s
@time con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 73 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)


