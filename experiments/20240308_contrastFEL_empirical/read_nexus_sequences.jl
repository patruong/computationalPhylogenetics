
function read_nexus(filename)
    # Open the file
    file = open(filename)

    # Initialize variables
    sequences = Dict()
    current_seq_name = ""

    # Read the file line by line
    for line in readline(file)
        # Skip empty lines and comments
        if isempty(line) || startswith(line, "#")
            continue
        end

        # Check if it's a data block
        if startswith(line, "BEGIN DATA;")
            continue
        end

        # Check if it's the end of the data block
        if startswith(line, "END;")
            break
        end

        # Check if it's a sequence name
        if startswith(line, ">")
            current_seq_name = strip(line, head=1)
            sequences[current_seq_name] = String("")
        else
            # Append the sequence data to the current sequence
            sequences[current_seq_name] *= line
        end
    end

    # Close the file
    close(file)

    # Return the dictionary of sequences
    return sequences
end

function read_file(file_path)
    open(file_path, "r") do file
        # Read the contents of the file into a string
        return read(file, String)
    end
end

function extract_sequences_nexus(nexus_content)
    seqnames = []
    seqs = []

    matrix_block = false
    if !isempty(nexus_content)
        for line in eachline(IOBuffer(nexus_content))
            # Skip empty lines and comments
            if isempty(line) || startswith(line, "#")
                continue
            end

            if matrix_block

                if line == "END;"
                    break
                end

                line = replace(line, r"\s+" => "")
                quoted_str = match(r"'([^']*)'", line).captures[1]
                quoted_str = replace(quoted_str, ";" => "")
                outside_str = replace(line, quoted_str => "")
                outside_str = replace(outside_str, "'" => "")
                outside_str = replace(outside_str, ";" => "")

                push!(seqnames, quoted_str)
                push!(seqs, outside_str)
            end

            if line == "MATRIX"
                matrix_block = true
            end

        end
    end
    return seqs, seqnames
end

function read_nexus_sequences(nexus_file)
    file_contents = read_file(nexus_file)
    seqs, seqnames = extract_sequences_nexus(file_contents)
    return seqs, seqnames
end

function write_seqarray_to_fasta(output, seq, seqnames)
    output_file = output
    open(output_file, "w") do f
        for i in 1:length(seqs)
            seq = seqs[i]
            seqname = seqnames[i]

            # Write sequence name
            write(f, ">$seqname\n")

            # Write sequence
            write(f, "$seq\n")
        end
    end
end

function extract_nexus_tree(nexus_content)
    tree_block = false
    if !isempty(nexus_content)
        for line in eachline(IOBuffer(nexus_content))
            # Skip empty lines and comments
            if isempty(line) || startswith(line, "#")
                continue
            end

            if tree_block

                if line == "END;"
                    break
                end
                line = replace(line, r"\s+" => "")
                line = replace(line, "TREEtree=" => "")
                return line
            end

            if line == "BEGIN TREES;"
                tree_block = true
            end

        end
    end

end

function read_nexus_tree(nexus_file)
    nexus_content = read_file(nexus_file)
    treestring = extract_nexus_tree(nexus_content)
    return treestring
end


function write_treestring(output, treestring)
    output_file = output
    open(output_file, "w") do f
        write(f, treestring)
    end
end

file_path = "/home/patrick/git/computationalPhylogenetics/data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex"
seqs, seqnames = read_nexus_sequences(file_path)
treestring = read_nexus_tree(file_path)

write_seqarray_to_fasta("/home/patrick/git/computationalPhylogenetics/experiments/20240308_contrastFEL_empirical/data/HIV_RT.fasta", seqs, seqnames)
write_treestring("/home/patrick/git/computationalPhylogenetics/experiments/20240308_contrastFEL_empirical/data/HIV_RT.nwk", treestring)

