using DataFrames
using CSV

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

	ID_ket = read_JBT("./exp/probe/ketamine/ketamine/") |> df_ID -> unique(df_ID.ID)

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
