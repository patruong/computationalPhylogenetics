# Define the parameters for the MG94 model
ω = 1.0  # Nonsynonymous/synonymous rate ratio
κ = 2.0  # Transition/transversion rate ratio
debug = true # debug message true or false

# Equilibrium codon frequencies for 61 sense codons (excluding stop codons)
# You need to specify the actual values for your specific dataset
# For example, you can calculate these frequencies from your data or use reference values.
#equilibrium_frequencies = [0.05, 0.05, ..., 0.05]  # Replace with your values

# Define the equilibrium frequencies for the 61 sense codons
equilibrium_frequencies = fill(1.0 / 61.0, 61)

# Create the rate matrix Q for the MG94 model
# The size of the matrix is 61x61 for the 61 sense codons
n_codons = length(equilibrium_frequencies)
Q = zeros(n_codons, n_codons)


function map_codon_to_base(codon::Int)
    non_stop_codons = [
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
        "TGG"
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
        "ATA" => "I", "ATC" => "I", "ATT" => "I", "ATG" => "M",
        "ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T",
        "AAC" => "N", "AAT" => "N", "AAA" => "K", "AAG" => "K",
        "AGC" => "S", "AGT" => "S", "AGA" => "R", "AGG" => "R",
        "CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L",
        "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P",
        "CAC" => "H", "CAT" => "H", "CAA" => "Q", "CAG" => "Q",
        "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R",
        "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V",
        "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A",
        "GAC" => "D", "GAT" => "D", "GAA" => "E", "GAG" => "E",
        "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G",
        "TCA" => "S", "TCC" => "S", "TCG" => "S", "TCT" => "S",
        "TTC" => "F", "TTT" => "F", "TTA" => "L", "TTG" => "L",
        "TAC" => "Y", "TAT" => "Y", "TAA" => "_", "TAG" => "_",
        "TGC" => "C", "TGT" => "C", "TGA" => "_", "TGG" => "W"
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
    pyrimidines = ['C', 'T']

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

            if debug == true
                println("i: $i")
                println("j: $j")
                println("codon_i: $codon_i")
                println("codon_j: $codon_j")
                println("count_difference: $(count_differences(codon_i, codon_j))")
                println("rate: $rate")
                println()
            end

            AA_i = map_codons_to_amino_acids([codon_i])
            AA_j = map_codons_to_amino_acids([codon_j])

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


# Fill the diagonal of the rate matrix with negative values
for i in 1:n_codons
    Q[i, i] = -sum(Q[i, :])
end


# The Q matrix now represents the rate matrix for the MG94 model

# To obtain equilibrium rates, you can exponentiate the Q matrix
t = 1.0  # Time (you can specify the desired time)
P_t = exp(Q * t)

# Equilibrium rates are given by the first row of the matrix P_t
equilibrium_rates = P_t[1, :]

# Print or use the equilibrium rates as needed
println("Equilibrium Rates: ", equilibrium_rates)

