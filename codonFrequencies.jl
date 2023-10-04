
module CodonFrequencies

export generate_codon_frequencies, generate_nucleotide_frequencies

function generate_nucleotide_frequencies(freq_list::Vector{Float64}, nucleotides::Vector{Char}=['A', 'C', 'G', 'T'])
    if length(freq_list) != length(nucleotides)
        throw(ArgumentError("Input lists must have the same length"))
    end

    nucleotide_frequencies = Dict{Char,Float64}()

    for (nucleotide, frequency) in zip(nucleotides, freq_list)
        nucleotide_frequencies[nucleotide] = frequency
    end

    return nucleotide_frequencies
end

function generate_codon_frequencies(freq_list::Vector{Float64}, nucleotide_list::Vector{Char}, remove_stop_codon::Bool=true, debug::Bool=false)
    nucleotide_frequencies = generate_nucleotide_frequencies(freq_list, nucleotide_list)

    # List of all codons (including stop codons)
    all_codons = [
        "ATA", "ATC", "ATT", "ATG",
        "ACA", "ACC", "ACG", "ACT",
        "AAC", "AAT", "AAA", "AAG",
        "AGC", "AGT", "AGA", "AGG",
        "CTA", "CTC", "CTG", "CTT",
        "CCA", "CCC", "CCG", "CCT",
        "CAC", "CAT", "CAA", "CAG",
        "CGA", "CGC", "CGG", "CGT",
        "GTA", "GTC", "GTG", "GTT",
        "GCA", "GCC", "GCG", "GCT",
        "GAC", "GAT", "GAA", "GAG",
        "GGA", "GGC", "GGG", "GGT",
        "TCA", "TCC", "TCG", "TCT",
        "TTC", "TTT", "TTA", "TTG",
        "TAC", "TAT", "TGC", "TGT",
        "TGG",
        "TAA", "TAG", "TGA"  # Stop codons
    ]

    # Calculate equilibrium frequencies for each codon
    codon_frequencies = Dict{String,Float64}()

    for codon in all_codons
        frequency = 1.0
        for nucleotide in codon
            frequency *= nucleotide_frequencies[nucleotide]
        end
        codon_frequencies[codon] = frequency
    end

    if debug == true
        # Print the equilibrium frequencies for each codon
        for (codon, frequency) in codon_frequencies
            println("Codon: $codon, Frequency: $frequency")
        end
    end

    if remove_stop_codon == true
        function if_not_stop_codon(codon)
            stop_codons = Set(["TAG", "TGA", "TAA", "UAA", "UGA", "UAG"])
            return !(codon in stop_codons)
        end

        filter!(kv -> if_not_stop_codon(kv[1]), codon_frequencies)
    end
    codon_frequencies = sort(collect(codon_frequencies))
    return codon_frequencies

end

function dict_to_list(dictionary)
    # Usage example 
    # Create a dictionary
    # my_dict = Dict(
    #     "A" => 1,
    #     "B" => 2,
    #     "C" => 3
    # )
    # keys_list, values_list = dict_to_list(my_dict)

    # Initialize empty lists for keys and values
    keys_list = []
    values_list = []

    # Extract keys and values from the dictionary
    for (key, value) in dictionary
        push!(keys_list, key)
        push!(values_list, value)
    end
    return keys_list, values_list
end

function generate_codon_frequencies_key_values(freq_list::Vector{Float64}, nucleotide_list::Vector{Char}, remove_stop_codon::Bool=true, debug::Bool=false)
    codon_frequency_dictionary = generate_codon_frequencies(freq_list, nucleotide_list, remove_stop_codon, debug)
    #println(codon_frequency_dictionary)
    codon_list, equilibrium_frequency_list = dict_to_list(codon_frequency_dictionary)
    return codon_list, equilibrium_frequency_list
end

end  # module



# Define nucleotide frequencies
# freq_list = [0.2, 0.3, 0.3, 0.2]
# nucleotide_list = ['A', 'C', 'G', 'T']
# debug = true
# codon_frequencies = generate_codon_frequencies(freq_list, nucleotide_list, debug)
# codon_frequencies


