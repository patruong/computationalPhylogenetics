using DataFrames

cd("/home/patrick/git/computationalPhylogenetics/experiments/20231231_results_check/difFUBAR_output_orig/")
file_baseline = "output/baseline/time_output.txt"
file_max = "output/max/time_output.txt"
file_patrick = "output/patrick/time_output.txt"
file_patrick_max = "output/patrick_max/time_output.txt"
file_patrick_max_child = "output/patrick_max_child/time_output.txt"
file_final = "output/final/time_output.txt"

function read_text_file(filename::AbstractString)
    file = open(filename, "r")
    content = read(file, String)
    close(file)
    return content
end

function parse_time_info(text::AbstractString)
    time_info = Dict{String,Float64}()
    pattern = r"(\w+)\s+(\d+m\d+\.\d+s)"
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
    pattern = r"(\w+)\s+(\d+m\d+\.\d+s)"
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

# Create an empty DataFrame
df = DataFrame(file=String[], real=Union{Float64,Missing}[], sys=Union{Float64,Missing}[], user=Union{Float64,Missing}[])

# Parse time for each file and add rows to the DataFrame
for file in [file_baseline, file_max, file_patrick, file_patrick_max, file_patrick_max_child, file_final]
    time_info = parse_time(file)
    push!(df, (file, time_info["real"], time_info["sys"], time_info["user"]))
end

# Display the resulting DataFrame
display(df)
