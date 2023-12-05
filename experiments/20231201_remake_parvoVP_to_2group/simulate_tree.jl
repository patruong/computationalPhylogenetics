using MolecularEvolution, Compose
using Cairo

cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231201_remake_parvoVP_to_2group")

tree = sim_tree(n=100)
#left
left_clade = getnodelist(tree.children[1])
#right
right_clade = getnodelist(tree.children[2])



function add_group_extension_str(clade, str)
    for node in clade
        node.name = (node.name * str)
    end
end

add_group_extension_str(left_clade, "{REFERENCE}")
add_group_extension_str(right_clade, "{TEST}")

img = tree_draw(tree, draw_labels=true, canvas_height=5cm)
img |> PDF("imgout_1.pdf", 10cm, 10cm)


newick_str = newick(tree)
file_path = "sim_tree.nwk"
file = open(file_path, "w")
write(file, newick_str)
close(file)


write_nexus("sim_tree.nex", tree)





