function loglikelihood_choices(model, chain)

	loglikelihoods = Turing.pointwise_loglikelihoods(model,
													MCMCChains.get_sections(chain,
																			:parameters))
	ynames = string.(keys(loglikelihoods))
	loglikelihoods_vals = getindex.(Ref(loglikelihoods), ynames)
	lr = permutedims(cat(loglikelihoods_vals...; dims=3), (2, 1, 3))
	return lr
end

function chain_to_idt(df, chain, func_R_history)

	(choices, data) = df_to_JBT(df, func_R_history)
	model = logistic_past(choices, data)

	loglikelihood_v = loglikelihood_choices(model, chain)

	idt = from_mcmcchains(chain;
					       log_likelihood=Dict("choices"=>loglikelihood_v),
					       library="Turing")
	return idt
end
