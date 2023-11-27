
for (row_ind, cp) in enumerate(codon_param_vec)
    alpha = cp[1]
    omegas = cp[2:end]
    tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

    #print("felsenstein")
    felsenstein!(tree, tagged_models)
    #This combine!() is needed because the current site_LLs function applies to a partition
    #And after a felsenstein pass, you don't have the eq freqs factored in.
    #We could make a version of log_likelihood() that returns the partitions instead of just the sum
    #print("combine")
    combine!.(tree.message, tree.parent_message)
    #println(length(tree.message[1]))
    #print("log_con_like_matrix")
    log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
    verbosity > 0 && if mod(row_ind, 500) == 1
        print(round(100 * row_ind / length(codon_param_vec)), "% ")
        flush(stdout)
    end
end