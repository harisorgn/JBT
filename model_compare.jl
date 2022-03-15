using Serialization
using LinearAlgebra
using ArviZ

include("JBT.jl")
include("read.jl")
include("metrics.jl")
include("model.jl")

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

#=
df = get_cohort_df(["amph"])
chain_learning = deserialize("chain_amph_learning.jls")
chain_RR = deserialize("chain_amph_RR.jls")
=#

df = get_batch_df(
				Dict(
					"./exp/probe/baseline/baseline_effect/"=>"baseline", 
					"./exp/probe/ketamine/ketamine_effect/"=>"ketamine", 
					"./exp/probe/ketamine/vehicle_effect/"=>"vehicle"
					)
				)

chain_RR = deserialize("chain_RR_effect_batch.jls")
chain_learning = deserialize("chain_learning_effect_batch.jls")

idt_learning = chain_to_idt(df, chain_learning, cue_reward_rate)
idt_RR = chain_to_idt(df, chain_RR, reward_rate)

d_model = Dict("learning"=>idt_learning, "RR"=>idt_RR)
compare(d_model, ic="loo", method="stacking", scale="log")
