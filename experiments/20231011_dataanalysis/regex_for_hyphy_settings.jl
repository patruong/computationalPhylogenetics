sim = 1 #1...530

file_name = "../../contrastFEL_data/omnibus/sims." * "$sim" * ".settings"
#json_data = JSON.parsefile(file_name)

file = open(file_name, "r")
str = readlines(file)
close(file)


i = 0
for line in str
    println(line)
    i += 1
    if i == 10
        break
    end
end


# Test the parser function on the given string
str = """{
 "model":{
   "EFV":    {
                {0.01929014233974128} 
                {0.01622765628659735} 
                {0.01769443322434847} 
                {0.01623235001696338} 
                {0.01622972861556908} 
                {0.01365311115696728} 
                {0.01488718145152558} 
                {0.0136570602190682} 
                {0.01769333756976861} 
                {0.01488435883925134} 
                {0.01622971850752687} 
                {0.01488866403069913} 
                {0.01623137334257151} 
                {0.01365449476855555} 
                {0.01488869012427688} 
                {0.01365844423085596} 
                {0.01929347190999687} 
                {0.0162304572572039}},
    "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
    "ID":"simulator.substitution_model",
    "Q":    [
        {"", "simulator.substitution_model.theta_AC*beta*0.23367778810339", "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", "simulator.substitution_model.theta_AT*beta*0.2337453776869083", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", ""}, 
        {"simulator.substitution_model.theta_AC*beta*0.2777775", "", "simulator.substitution_model.theta_CG*beta*0.2547993342097017", "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", ""}, 
        {"simulator.substitution_model.theta_AG*alpha*0.2777775", "simulator.substitution_model.theta_CG*beta*0.23367778810339", ""} 
        ]
 }
}"""

# Lets try regex clean up this data.

str_ = replace(str, " " => "")
str_ = replace(str_, "\n" => "")
# Use regex to replace `{` with `[` when not in the format '{"some text here":'
str_ = replace(str_, r"(?<!{\"[^}]+\":)\{" => "[")

println(str_)


# Use regex to replace ":[" with ":[["
str = replace(str, r"(\":\[\s*)" => s -> s * "[")
str = replace(str, r"(\":\[\s*)" => s -> s * "[")
str = replace(str, r"\],\}" => "]],}")
str = replace(str, r"\]\],\}" => "]]}")

