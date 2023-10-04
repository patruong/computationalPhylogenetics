include("codonFrequencies.jl")  # Include the script to load the module

import .CodonFrequencies.generate_codon_frequencies_key_values
import .CodonFrequencies.generate_codon_frequencies

module MG94


function map_codon_to_base(codon::Int)
    non_stop_codons = [
        "AUA", "AUC", "AUU", "AUG",
        "ACA", "ACC", "ACG", "ACU",
        "AAC", "AAU", "AAA", "AAG",
        "AGC", "AGU", "AGA", "AGG",
        "CUA", "CUC", "CUG", "CUU",
        "CCA", "CCC", "CCG", "CCU",
        "CAC", "CAU", "CAA", "CAG",
        "CGA", "CGC", "CGG", "CGU",
        "GUA", "GUC", "GUG", "GUU",
        "GCA", "GCC", "GCG", "GCU",
        "GAC", "GAU", "GAA", "GAG",
        "GGA", "GGC", "GGG", "GGU",
        "UCA", "UCC", "UCG", "UCU",
        "UUC", "UUU", "UUA", "UUG",
        "UAC", "UAU", "UGC", "UGU",
        "UGG"
    ]

    if codon < 1 || codon > 61
        throw(ArgumentError("Codon must be in the range 1 to 61"))
    end

    return non_stop_codons[codon]
end

function map_codons_to_amino_acids(codon_list)
    # Define a dictionary to map codons to amino acids
    # Example usage:
    # non_stop_codons = [
    #     "ATG", "ATT", "ATT", "ATA", "ATG", "CTT", "CTT", "CTA", "CTG", "ATT",
    #     "ATC", "ATA", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG",
    #     "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC",
    #     "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC",
    #     "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT",
    #     "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA",
    #     "GGG"
    # ]
    # 
    # amino_acids = map_codons_to_amino_acids(non_stop_codons)
    # println(amino_acids)
    codon_to_amino_acid = Dict(
        "AUA" => "I", "AUC" => "I", "AUU" => "I", "AUG" => "M",
        "ACA" => "U", "ACC" => "T", "ACG" => "U", "ACU" => "U",
        "AAC" => "N", "AAU" => "N", "AAA" => "K", "AAG" => "K",
        "AGC" => "S", "AGU" => "S", "AGA" => "R", "AGG" => "R",
        "CUA" => "L", "CUC" => "L", "CUG" => "L", "CUU" => "L",
        "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCU" => "P",
        "CAC" => "H", "CAU" => "H", "CAA" => "Q", "CAG" => "Q",
        "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGU" => "R",
        "GUA" => "V", "GUC" => "V", "GUG" => "V", "GUU" => "V",
        "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCU" => "A",
        "GAC" => "D", "GAU" => "D", "GAA" => "E", "GAG" => "E",
        "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGU" => "G",
        "UCA" => "S", "UCC" => "S", "UCG" => "S", "UCU" => "S",
        "UUC" => "F", "UUU" => "F", "UUA" => "L", "UUG" => "L",
        "UAC" => "Y", "UAU" => "Y", "UAA" => "_", "UAG" => "_",
        "UGC" => "C", "UGU" => "C", "UGA" => "_", "UGG" => "W"
    )

    # Initialize an empty array to store the amino acids
    amino_acids = []

    # Iterate through the codon list and map each codon to its amino acid
    for codon in codon_list
        amino_acid = get(codon_to_amino_acid, codon, "Unknown")
        push!(amino_acids, amino_acid)
    end

    return amino_acids
end

function count_differences(seq1::String, seq2::String)
    if length(seq1) != length(seq2)
        throw(ArgumentError("Strings must have the same length"))
    end

    count = 0
    for i in 1:length(seq1)
        if seq1[i] != seq2[i]
            count += 1
        end
    end

    return count
end


function is_transition(base1::Char, base2::Char)
    # Define purines and pyrimidines
    purines = ['A', 'G']
    pyrimidines = ['C', 'T', 'U']

    # Check if the bases are valid
    if !(base1 in purines ∪ pyrimidines) || !(base2 in purines ∪ pyrimidines)
        throw(ArgumentError("Invalid bases"))
    end

    if (base1 in purines && base2 in purines) || (base1 in pyrimidines && base2 in pyrimidines)
        return true
    else
        return false
    end
end

function is_transversion(base1::Char, base2::Char)
    # A transversion is any nucleotide substitution that is not a transition
    return !is_transition(base1, base2)
end

function codon_position_diff(codon_i, codon_j)
    # Convert strings to character arrays
    c1 = collect(codon_i)
    c2 = collect(codon_j)

    # Find the first differing position using broadcasting
    position = findfirst(c1 .!= c2)
end

function compute_Q(equilibrium_frequencies::Vector, α::Float64, β::Float64, κ::Float64, debug::Bool=false)
    # Q are the instantaneous rates
    # Create the rate matrix Q for the MG94 model
    # The size of the matrix is 61x61 for the 61 sense codons
    n_codons = length(equilibrium_frequencies)
    Q = zeros(n_codons, n_codons)
    ω = β / α

    # Populate the rate matrix Q based on the MG94 model parameters
    for i in 1:n_codons
        for j in 1:n_codons
            if i != j
                codon_i = map_codon_to_base(i)
                codon_j = map_codon_to_base(j)
                # Calculate the instantaneous rate for transition or transversion
                rate = 1.0  # Default rate (for non-identical codons)
                # If there is more than one substitution set rate to 0
                if count_differences(codon_i, codon_j) > 1
                    rate = 0
                else # else check at which position the difference is 
                    pos = codon_position_diff(codon_i, codon_j)
                    if is_transition(codon_i[pos], codon_j[pos])
                        rate = κ
                    elseif is_transversion(codon_i[pos], codon_j[pos])
                        rate = 1.0
                    end
                end

                AA_i = map_codons_to_amino_acids([codon_i])
                AA_j = map_codons_to_amino_acids([codon_j])

                if debug == true
                    println("i: $i")
                    println("j: $j")
                    println("codon_i: $codon_i")
                    println("codon_j: $codon_j")
                    println("count_difference: $(count_differences(codon_i, codon_j))")
                    println("AA_i: $(AA_i[1])")
                    println("AA_j: $(AA_j[1])")
                    println("rate: $rate")
                    println()
                end

                # Populate the rate matrix Q
                if AA_i == AA_j
                    synonymous_rate = rate * equilibrium_frequencies[j]
                    Q[i, j] = synonymous_rate
                else
                    nonsynonymous_rate = rate * equilibrium_frequencies[j] * ω
                    Q[i, j] = nonsynonymous_rate
                end
            end
        end
    end


    # Q_ij are the instantaneous rates
    # Fill the diagonal of the rate matrix with negative values
    for i in 1:n_codons
        Q[i, i] = -sum(Q[i, :])
    end

    return Q
end

function compute_equilibrium_rate(Q::Matrix{Float64}, t::Float64=1.0)
    # To obtain equilibrium rates, you can exponentiate the Q matrix
    P_t = exp(Q * t)

    # Equilibrium rates are given by the first row of the matrix P_t
    equilibrium_rates = P_t[1, :]
    return equilibrium_rates
end


function compute_Q(equilibrium_frequencies::Vector, α::Float64, β::Float64, κ::Float64, debug::Bool=false)
    # Q are the instantaneous rates
    # Create the rate matrix Q for the MG94 model
    # The size of the matrix is 61x61 for the 61 sense codons
    n_codons = length(equilibrium_frequencies)
    Q = zeros(n_codons, n_codons)
    ω = β / α

    # Populate the rate matrix Q based on the MG94 model parameters
    for i in 1:n_codons
        for j in 1:n_codons
            if i != j
                codon_i = map_codon_to_base(i)
                codon_j = map_codon_to_base(j)
                # Calculate the instantaneous rate for transition or transversion
                rate = 1.0  # Default rate (for non-identical codons)
                # If there is more than one substitution set rate to 0
                if count_differences(codon_i, codon_j) > 1
                    rate = 0
                else # else check at which position the difference is 
                    pos = codon_position_diff(codon_i, codon_j)
                    if is_transition(codon_i[pos], codon_j[pos])
                        rate = κ
                    elseif is_transversion(codon_i[pos], codon_j[pos])
                        rate = 1.0
                    end
                end

                AA_i = map_codons_to_amino_acids([codon_i])
                AA_j = map_codons_to_amino_acids([codon_j])

                if debug == true
                    println("i: $i")
                    println("j: $j")
                    println("codon_i: $codon_i")
                    println("codon_j: $codon_j")
                    println("count_difference: $(count_differences(codon_i, codon_j))")
                    println("AA_i: $(AA_i[1])")
                    println("AA_j: $(AA_j[1])")
                    println("rate: $rate")
                    println()
                end

                # Populate the rate matrix Q
                if AA_i == AA_j
                    synonymous_rate = rate * equilibrium_frequencies[j]
                    Q[i, j] = synonymous_rate
                else
                    nonsynonymous_rate = rate * equilibrium_frequencies[j] * ω
                    Q[i, j] = nonsynonymous_rate
                end
            end
        end
    end


    # Q_ij are the instantaneous rates
    # Fill the diagonal of the rate matrix with negative values
    for i in 1:n_codons
        Q[i, i] = -sum(Q[i, :])
    end

    return Q
end







end  # module


# # Define the parameters for the MG94 model
# α = 1.0 # synonymous rate 
# β = 1.0 # Nonsynonymous rate
# #ω = 1.0  # Nonsynonymous/synonymous rate ratio
# ω = β / α
# κ = 2.0  # Transition/transversion rate ratio
# debug = true # debug message true or false
# remove_stop_codon = true

# freq_list = [0.2, 0.3, 0.3, 0.2] # Nucleotide frequencies 
# nucleotide_list = ['A', 'C', 'G', 'T'] # Nucleotides
# codon_frequency_dictionary = generate_codon_frequencies(freq_list, nucleotide_list, remove_stop_codon, debug)


# # Equilibrium codon frequencies for 61 sense codons (excluding stop codons)
# # You need to specify the actual values for your specific dataset
# # For example, you can calculate these frequencies from your data or use reference values.
# #equilibrium_frequencies = [0.05, 0.05, ..., 0.05]  # Replace with your values

# # Define the equilibrium frequencies for the 61 sense codons
# equilibrium_frequencies = fill(1.0 / 61.0, 61)

# # These list includes stop frequencies
# codon_list, codon_freq_list = generate_codon_frequencies_key_values(freq_list, nucleotide_list, remove_stop_codon, debug)

# equilibrium_frequencies = codon_freq_list



# # The Q matrix now represents the rate matrix for the MG94 model
# Q = compute_Q(equilibrium_frequencies, α, β, κ, false)

# t = 1.0  # Time (you can specify the desired time)
# equilibrium_rates = compute_equilibrium_rate(Q, t)
# # Print or use the equilibrium rates as needed
# println("Equilibrium Rates: ", equilibrium_rates)

