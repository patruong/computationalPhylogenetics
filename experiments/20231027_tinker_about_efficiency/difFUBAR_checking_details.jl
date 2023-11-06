using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution


cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
analysis_name = "nobackground/Ace2"
seqnames, seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=false, verbosity=1)
@time code = MolecularEvolution.universal_code

seqnames
seqs

@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s

tree
GTRmat
F3x4_freqs
eq_freqs

#difFUBAR grid
con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s


function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end


tr(x) = 10^x - 0.05
trinv(x) = log10(x + 0.05)
#Defines the grid used for inference.
function gridsetup(lb, ub, num_below_one, trin, tr)
    step = (trin(1.0) - trin(lb)) / num_below_one
    return tr.(trin(lb):step:trin(ub))
end

foreground_grid = 6
background_grid = 4
alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2
num_groups = length(tags)
[println(n.name) for n in getnodelist(tree)]
[tags for n in getnodelist(tree)]
is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
tensor_dims = 1 + num_groups + is_background

function add_to_each_element(vec_of_vec, elems)
    return [vcat(v, [e]) for v in vec_of_vec for e in elems]
end

# check what is happening in the loop
codon_param_vec
omegagrid
a = [[1], [2], [3]]
b = [1, 2, 3]
a = add_to_each_element(a, b)
a = add_to_each_element(a, b)


codon_param_vec = [[a] for a in alphagrid]
param_kinds = ["Alpha"]

for g in 1:num_groups
    push!(param_kinds, "OmegaG$(g)")
    codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
end

param_kinds
codon_param_vec


@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=false)



@assert size(df) == (19, 8)


alpha
beta

# rebuilding felsenstein
using Plots, Phylo, Compose

cached_model = CodonMolecularEvolution.MG94_cacher(code)


function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
    models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
    return n::FelNode -> [models[model_ind(n.name, tags)]]
end

omegas = codon_param_vec[1]


tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)


felsenstein!(tree, tagged_models)


d = marginal_state_dict(tree, tagged_models);

compose_dict = Dict()
for n in getnodelist(tree)
    compose_dict[n] = (x, y) -> pie_chart(x, y, d[n][1].state[:, 1], size=0.02, opacity=0.75)
end
img = tree_draw(tree, draw_labels=false, line_width=0.5mm, compose_dict=compose_dict)


### draw test
str_tree = "(((((tax24:0.09731668728575642,(tax22:0.08792233964843627,tax18:0.9210388482867483):0.3200367900275155):0.6948314526087965,(tax13:1.9977212308725611,(tax15:0.4290074347886068,(tax17:0.32928401808187824,(tax12:0.3860215462534818,tax16:0.2197134841232339):0.1399122681886174):0.05744611946245004):1.4686085778061146):0.20724159879522402):0.4539334554156126,tax28:0.4885576926440158):0.002162260013924424,tax26:0.9451873777301325):3.8695419798779387,((tax29:0.10062813251515536,tax27:0.27653633028085006):0.04262434258357507,(tax25:0.009345653929737636,((tax23:0.015832941547076644,(tax20:0.5550597590956172,((tax8:0.6649025646927402,tax9:0.358506423199849):0.1439516404012261,tax11:0.01995439013213013):1.155181296134081):0.17930021667907567):0.10906638146207207,((((((tax6:0.013708993438720255,tax5:0.061144001556547097):0.1395453591567641,tax3:0.4713722705245479):0.07432598428904214,tax1:0.5993347898257291):1.0588025698844894,(tax10:0.13109032492533992,(tax4:0.8517302241963356,(tax2:0.8481963081549965,tax7:0.23754095940676642):0.2394313086297733):0.43596704123297675):0.08774657269409454):0.9345533723114966,(tax14:0.7089558245245173,tax19:0.444897137240675):0.08657675809803095):0.01632062723968511,tax21:0.029535281963725537):0.49502691718938285):0.25829576024240986):0.7339777396780424):4.148878039524972):0.0"
newt = gettreefromnewick(str_tree, FelNode)
ladderize!(newt)
compose_dict = Dict()
for n in getleaflist(newt)
    #Replace the rand(4) with the frequencies you actually want.
    compose_dict[n] = (x, y) -> pie_chart(x, y, MolecularEvolution.sum2one(rand(4)), size=0.03)
end
tree_draw(newt, draw_labels=false, line_width=0.5mm, compose_dict=compose_dict)


img = tree_draw(newt)
display(img)
img |> SVG("imgout.svg", 10cm, 10cm)
using Cairo
img |> PDF("imgout.pdf", 10cm, 10cm)



