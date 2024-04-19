using Glob
cd("/home/patrick/git/computationalPhylogenetics/contrastFEL_data/omnibus")

function read_fasta(filename::String)
    fasta_dict = Dict{String,String}()
    current_id = ""
    current_sequence = ""

    for line in eachline(filename)
        if startswith(line, '>')
            # If a line starts with '>', it indicates a new sequence
            if !isempty(current_id)
                fasta_dict[current_id] = current_sequence
            end
            current_id = split(line[2:end])[1]  # Extracting the ID
            current_sequence = ""
        else
            current_sequence *= strip(line)  # Concatenate sequence lines
        end
    end

    # Adding the last sequence to the dictionary
    if !isempty(current_id)
        fasta_dict[current_id] = current_sequence
    end

    return fasta_dict
end

function sequence_count(sequences)
    all_sequences = values(sequences)
    num_all_seq = length(all_sequences)
    unique_sequences = unique(all_sequences)
    num_unique_sequences = length(unique_sequences)
    ratio = num_unique_sequences / num_all_seq
    #println("Number of all sequences:", num_all_seq)
    #println("Number of unique sequences:", num_unique_sequences)
    #println("Ratio of unique seqeunces to all sequences:", ratio)
    return num_all_seq, num_unique_sequences, ratio
end

files = readdir()
filtered_files = filter(file -> occursin(r"[1-5]$", file), files)
file = filtered_files[1]
fasta_file = file
parts = split(fasta_file, '.')
sim = parts[2]
rep = parts[5]
sequences = read_fasta(fasta_file)
n_all_seq, n_unique_seq, all_to_unique_ratio = sequence_count(sequences)


