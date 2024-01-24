using DataFrames
df = DataFrame(a=1:7, b=4:10, c=7:13)
df[!, "test test"] = df[!, :a]
df[!, "test"] = df[!, :a]

function calculate_result(x, y, df)
    a_sum = sum(df[1:x, :a])
    return a_sum * y
end

function test2(x, y)
    return 2x + y
end
transform(df, AsTable([:a, :c]) => ByRow(x -> (calculate_result(x.a, x.c, df), 4x.c)) => ["ae what", "f"])

