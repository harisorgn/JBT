using Serialization

include("JBT.jl")
include("read.jl")
include("metrics.jl")
include("model.jl")

function run_cohort_ketamine(func_R_history, model_label)

	df = get_cohort_df(["1vs1", "4vs1", "ket"], 
						ID_excluded=["HO2_4", "HO2_11"], 
						S_excluded=["1vs1_1"])

 	(choices, data) = df_to_JBT(df, func_R_history)
	model = logistic_past(choices, data)
	chain = sample(model, NUTS(), MCMCThreads(), 2000, 4)

	serialize(string("chain_ket_", model_label, ".jls"), chain)
end

function run_cohort_amphetamine(func_R_history, model_label)

	df = get_cohort_df(["amph"])
	(choices, data) = df_to_JBT(df, func_R_history)
	model = logistic_past(choices, data)
	chain = sample(model, NUTS(), MCMCThreads(), 2000, 4)

	serialize(string("chain_amph_", model_label, ".jls"), chain)
end

function run_batch_ketamine_no_effect(func_R_history, model_label)

	d = Dict(
			"./exp/probe/baseline/baseline_no_effect/"=>"baseline", 
			"./exp/probe/ketamine/ketamine_no_effect/"=>"ketamine", 
			"./exp/probe/ketamine/vehicle_no_effect/"=>"vehicle"
			)

	df = get_batch_df(d)
	(choices, data) = df_to_JBT(df, func_R_history)
	model = logistic_past(choices, data)
	chain = sample(model, NUTS(), MCMCThreads(), 1000, 4)

	serialize(string("chain_", model_label, "_no_effect_batch.jls"), chain)
end

function run_batch_ketamine_effect(func_R_history, model_label)

	d = Dict(
			"./exp/probe/baseline/baseline_effect/"=>"baseline", 
			"./exp/probe/ketamine/ketamine_effect/"=>"ketamine", 
			"./exp/probe/ketamine/vehicle_effect/"=>"vehicle"
			)

	df = get_batch_df(d)
	(choices, data) = df_to_JBT(df, func_R_history)
	model = logistic_past(choices, data)
	chain = sample(model, NUTS(), MCMCThreads(), 1000, 4)

	serialize(string("chain_", model_label, "_effect_batch.jls"), chain)
end

Turing.setadbackend(:reversediff)

run_cohort_amphetamine(reward_rate, "RR")
run_cohort_amphetamine(cue_reward_rate, "learning")
