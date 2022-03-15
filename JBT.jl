
struct JBT
	ID::Vector{Int64}
	conditions::Vector{Int64}
	cues::Vector{Float64}
	R::Vector{Float64}
	I_m1::Vector{Float64}
	n_subjects::Int64
	n_conditions::Int64
end

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

	return (choices, JBT(ID, conditions, cues, R, I_m1, length(unique(df.ID)), length(unique(conditions))))
end

function get_cohort_df(study_v; ID_excluded=[], S_excluded=[])

	df = DataFrame(
					ID = String[], 
					trial = Int64[],
					date = String[],
					session = String[], 
					cue = Int64[], 
					choice = Int64[], 
					RT = Float64[], 
					outcome = Int64[]
					)

	for study in study_v

		dir = string("./exp/HO2/", study,"/")
		df_study = read_JBT(dir)

		filter!(row -> row.session == "probe", df_study)

		df_cb = CSV.File(string(dir,"counterbalance/", study, "_counterbalance.csv")) |> DataFrame
		df_study.session = map(row -> df_cb[df_cb.ID.==row.ID, row.date][1], eachrow(df_study))

		append!(df, df_study)
	end

	filter!(row -> (row.choice!=0) && 
					(row.cue!=0) && 
					!in(row.ID, ID_excluded) && 
					!in(row.session, S_excluded), 
			df)
	return df
end

function get_batch_df(study_dict)

	ID_ket = read_JBT("./exp/probe/ketamine/ketamine/") |> df_ID->unique(df_ID.ID)

	df = DataFrame(
					ID = String[], 
					trial = Int64[],
					date = String[],
					session = String[], 
					cue = Int64[], 
					choice = Int64[], 
					RT = Float64[], 
					outcome = Int64[]
					)

	for (dir, condition_name) in study_dict

		df_study = read_JBT(dir)

		filter!(row -> row.session == "probe", df_study)

		df_study.session = fill(condition_name, nrow(df_study))

		append!(df, df_study)
	end

	filter!(row -> (row.choice!=0) && (row.cue!=0) && in(row.ID, ID_ket), df)
	return df
end

condition_dict(df) = Dict(k=>i for (k,i) in zip(unique(df.session),eachindex(unique(df.session))))

