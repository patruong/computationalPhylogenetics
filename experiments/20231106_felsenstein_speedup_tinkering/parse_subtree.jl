
# A function that takes a Newick string and returns a list of sub-trees
function extract_subtrees(newick)
    # Initialize an empty list to store the sub-trees
    subtrees = []
    # Initialize a stack to keep track of the parentheses
    stack = []
    # Initialize a variable to store the current group
    group = ""
    # Initialize a variable to store the current sub-tree
    subtree = ""
    # Initialize a variable to store the current index in the Newick string
    i = 1
    # Loop through each character in the Newick string
    while i <= length(newick)
        c = newick[i]
        # If the character is an opening parenthesis, push it to the stack
        if c == '('
            push!(stack, c)
            # If the character is a closing parenthesis, pop from the stack
        elseif c == ')'
            pop!(stack)
            # If the character is a curly bracket, extract the group name
        elseif c == '{'
            # Find the index of the matching closing bracket
            j = findnext('}', newick, i)
            # Get the group name between the brackets
            group = newick[i+1:j-1]
            # Skip the rest of the characters until the closing bracket
            i = j
            continue
            # If the character is a semicolon, end the loop
        elseif c == ';'
            break
        end
        # Append the character to the current sub-tree
        subtree *= c
        # If the stack is empty, it means we have reached the end of a sub-tree
        if isempty(stack)
            # If the group name is not empty, it means we have a valid sub-tree
            if group != ""
                # Append the group name and a semicolon to the sub-tree
                subtree *= "{$group};"
                # Add the sub-tree to the list of sub-trees
                push!(subtrees, subtree)
                # Reset the group name and the sub-tree
                group = ""
                subtree = ""
            end
        end
        # Increment the index by one
        i += 1
    end
    # Return the list of sub-trees
    return subtrees
end

# A sample input
newick = "((((N156{TEST}:0.05022147571229687,N157{TEST}:0.05641084875675355)N293{TEST}:0.1737237766768709,(N158{TEST}:0.21803248337425,N159{TEST}:0.1194040491539744)N294{REFERENCE}:0.3313542119310432)N361{TEST}:0.08722508335508614,((N160{TEST}:0.002646079914739969,N161{TEST}:0.05502129376593578)N295{REFERENCE}:0.1035036408828642,(N162{REFERENCE}:0.01210067424973002,N163{TEST}:0.1221041986807102)N296{TEST}:0.02891981622318314)N362{TEST}:0.1591913049511287)N395{TEST}:0.06117152083537829,(((N164{TEST}:0.05457434945527553,N165{TEST}:0.03411432503522936)N297{REFERENCE}:0.06628947303271021,(N166{TEST}:0.08681644831599826,N167{TEST}:0.2031348023473614)N298{TEST}:0.03127432538177079)N363{TEST}:0.01794812319423736,((N168{TEST}:0.2021384493771804,N169{TEST}:0.1422922322989518)N299{TEST}:0.05838849921499194,(N170{TEST}:0.07697879124598514,N171{TEST}:0.06243846792898091)N300{TEST}:0.2339237387089033)N364{TEST}:0.01397715112067482)N396{TEST}:0.2065996151697533):0.0;"


pattern = r"\((.*?)\)"
str = "This is a (test) string (results)."

match(pattern, str)
match_vector = findall(pattern, str)

str[match_vector[1]]
str[match_vector[2]]

i = 1
while i <= length(newick)
    c = newick[i]
    println(c)
    i += 1
end




# Call the function and print the output
#output = extract_subtrees(input)
#for subtree in output
#    println(subtree)
#end
#
#output
#
#show(output)


