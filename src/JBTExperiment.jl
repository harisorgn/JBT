struct JBTExperiment
	ID::Vector{Int64}
	conditions::Vector{Int64}
	cues::Vector{Float64}
	R::Vector{Float64}
	I_m1::Vector{Float64}
	n_subjects::Int64
	n_conditions::Int64
end

condition_dict(df) = Dict(k=>i for (k,i) in zip(unique(df.session),eachindex(unique(df.session))))

function df_to_JBT(df, func_R_history)

	d_ID = Dict(k=>i for (k,i) in zip(unique(df.ID),eachindex(unique(df.ID))))
	d_condition = Dict(k=>i for (k,i) in zip(unique(df.session),eachindex(unique(df.session))))

	choices = Int64[]
	conditions = Int64[]
	ID = Int64[]
	cues = Float64[]
	R = Float64[]
	I_m1 = Float64[]

	for df_condition in groupby(df, :session)
		for df_ID in groupby(df_condition, :ID)

			@assert length(unique(df_ID.session))==1

			mask_amb = append!([false], df_ID.cue[2:end].==5)

			append!(choices, map(c -> c==8 ? 1 : 0, df_ID.choice[mask_amb]))

			n_amb_trials = count(x -> x==5, df_ID.cue[mask_amb])

			append!(cues, fill(5, n_amb_trials))
			append!(ID, fill(d_ID[df_ID.ID[1]], n_amb_trials))
			append!(conditions, fill(d_condition[df_ID.session[1]], n_amb_trials))

			df_R = transform(df_ID, AsTable([:cue, :outcome, :trial]) => (x -> func_R_history(x)) => :R)
			append!(R, df_R.R[mask_amb])

			df_I_m1 = transform(df_ID, AsTable([:choice, :outcome, :session]) => (x -> wsls_trial_index(x)) => :I)
			append!(I_m1, df_I_m1.I[mask_amb])
		end
	end

	return (choices, JBTExperiment(ID, conditions, cues, R, I_m1, length(unique(df.ID)), length(unique(conditions))))
end
