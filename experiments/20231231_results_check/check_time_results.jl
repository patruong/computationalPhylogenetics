using DataFrames
using CSV

cd("/home/patrick/git/computationalPhylogenetics/experiments/20231231_results_check/difFUBAR_output_orig/")
file_baseline = "output/baseline/time_output.txt"
file_max = "output/max/time_output.txt"
file_patrick = "output/patrick/time_output.txt"
file_patrick_max = "output/patrick_max/time_output.txt"
file_patrick_max_child = "output/patrick_max_child/time_output.txt"
file_final = "output/final/time_output.txt"

function read_text_file(filename::AbstractString)
    # Open the file in read mode
    file = open(filename, "r")

    # Read the contents of the file
    content = read(file, String)

    # Close the file
    close(file)

    return content
end

function parse_time_info(text::AbstractString)
    time_info = Dict{String,Float64}()

    # Regular expression pattern to match lines with time information
    pattern = r"(\w+)\s+(\d+m\d+\.\d+s)"

    # Iterate through lines and extract time information
    for line in eachline(IOBuffer(text))
        match_result = match(pattern, line)
        if match_result !== nothing
            key, value = match_result.captures
            minutes, seconds = match(r"(\d+)m(\d+\.\d+)s", value).captures
            total_seconds = parse(Float64, minutes) * 60 + parse(Float64, seconds)
            time_info[key] = total_seconds
        end
    end

    return time_info
end

function parse_time_info_safe(text::AbstractString)
    time_info = Dict{String,Union{Float64,Missing}}()

    # Regular expression pattern to match lines with time information
    pattern = r"(\w+)\s+(\d+m\d+\.\d+s)"

    # Iterate through lines and extract time information
    for line in eachline(IOBuffer(text))
        try
            match_result = match(pattern, line)
            if match_result !== nothing
                key, value = match_result.captures
                minutes, seconds = match(r"(\d+)m(\d+\.\d+)s", value).captures
                total_seconds = parse(Float64, minutes) * 60 + parse(Float64, seconds)
                time_info[key] = total_seconds
            end
        catch
            # Set missing values in case of an error
            time_info[key] = missing
        end
    end

    return time_info
end


function parse_time(file)
    try
        return parse_time_info_safe(read_text_file(file))
    catch
        time_info = Dict{String,Union{Float64,Missing}}()
        time_info["real"] = missing
        time_info["sys"] = missing
        time_info["user"] = missing
        return time_info
    end
end

parse_time(file_baseline)
parse_time(file_max)
parse_time(file_patrick)
parse_time(file_patrick_max)
parse_time(file_patrick_max_child)
parse_time(file_final)

