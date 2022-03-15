using DataFrames
using CSV

abstract type session_t end
struct probe_t <: session_t end
struct probe_1vs1_t <: session_t end
struct probe_compound_1vs1_t <: session_t end
struct train_1vs1_t <: session_t end
struct train_4vs1_t <: session_t end
struct train_1vs1_compound_t <: session_t end

function date_sort!(fv::Array{String})

	month_d = Dict(
					"Jan"=>1,
					"Feb"=>2,
					"Mar"=>3,
					"Apr"=>4,
					"May"=>5,
					"Jun"=>6,
					"Jul"=>7,
					"Aug"=>8,
					"Sep"=>9,
					"Oct"=>10,
					"Nov"=>11,
					"Dec"=>12
					)

	date_v = map(fv) do date
		(day, month, year) = split(date,"-")
		day = tryparse(Int64, day)
		month = month_d[month]
		year = tryparse(Int64, year[1:4])
		(day, month, year, date)
	end

	sort!(date_v, by = x->(x[3],x[2],x[1],x[4]))
	fv[:] = map(x -> x[4], date_v) 
end

function read_csv_var_cols(file_path::String)

	# read a csv file with a variable number of columns across rows
	# missing D elements are filled with ""

	max_ncols = 0 
	nrows = 0 

	open(file_path) do f
		while !eof(f)
			new_line = split(strip(readline(f)),',') 
			if length(new_line) > max_ncols
				max_ncols = length(new_line)
			end
			nrows += 1
		end
	end

	D = Matrix{String}(undef, nrows, max_ncols)

	open(file_path) do f

		headers = [string("col",i) for i = 1 : max_ncols]
		ncols = max_ncols

		i = 1
		while !eof(f)

			new_line = split(strip(readline(f)),',')
			length(new_line) < ncols && append!(new_line,["" for i=1:ncols-length(new_line)])

			D[i, :] = new_line

			i += 1
		end
	end
	return D
end

function read_JBT(dir; reversed=false)

	stype_d = Dict(
					"Probe tone midpoint only testing" => (probe_t(), 20),
					"1 vs 1 Probe Test" => (probe_1vs1_t(), 22),
					"Discrimination training" => (train_1vs1_t(), 18),
					"Pure tones testing" => (train_4vs1_t(), 18),
					"1vs1" => (probe_compound_1vs1_t(), 32),
					"Discrimination" => (train_1vs1_compound_t(), 42)
					)

	file_v = filter(x->occursin(".csv", x), readdir(dir))
	date_sort!(file_v)

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

	for (file_name, session) in zip(file_v, eachindex(file_v))

		D = read_csv_var_cols(string(dir, file_name)) 

		(nrows, ~) = size(D)	

		ID_v = D[occursin.("Id",D[:,1]), 2]

		s = D[occursin.("AC Comment",D[:,1]), 3]
		stype_v = map(x -> stype_d[x][1], s)	
		ncols_v = map(x -> stype_d[x][2], s)

		idx_staRT_v = findall(x->occursin("Outcome",x), D[:,2]) .+ 1
		idx_end_v = findall(x->occursin("ENDDATA",x), D[:,1]) .- 1

		@assert length(idx_staRT_v) == length(idx_end_v)

		for (idx_s, idx_e, i) in zip(idx_staRT_v, idx_end_v, eachindex(idx_staRT_v))

			subj_D = parse.(Int64, D[idx_s:idx_e, 1:ncols_v[i]])

			append_subj!(df, subj_D, ID_v[i], split(file_name,"_")[1], stype_v[i]; 
						reversed=reversed)
		end
	end		

	return df
end

function append_subj!(
					df, 
					D, 
					subj_id::String, 
					date, 
					::probe_t; 
					col_offset=0,
					reversed=false
					)

	RT_max = 20.0 # sec
	n_trials = size(D,1)

	outcome = zeros(Int64, n_trials)
	choice = zeros(Int64, n_trials)
	cue = zeros(Int64, n_trials)
	RT = zeros(Float64, n_trials)
	
	id_number_idx = findlast(isequal('_'), subj_id)
	id_number = tryparse(Int64, subj_id[id_number_idx+1:end])

	# p1 : probe cue playing and route 1 lever (left) was set as correct
	# p2 : probe cue playing and route 2 lever (right) was set as correct

	mask_corr_cue1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,13 + col_offset])
	mask_corr_cue2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,15 + col_offset])
	mask_incorr_cue1 = map((x,y) -> x == 1 && y != 0, D[:,2], D[:,13 + col_offset])
	mask_incorr_cue2 = map((x,y) -> x == 3 && y != 0, D[:,2], D[:,15 + col_offset])
	mask_corr_p1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,17 + col_offset])
	mask_corr_p2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,19 + col_offset])
	mask_incorr_p1 = map((x,y) -> x == 1 && y != 0, D[:,2], D[:,17 + col_offset])
	mask_incorr_p2 = map((x,y) -> x == 3 && y != 0, D[:,2], D[:,19 + col_offset])
	mask_om_cue1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 13 + col_offset])
	mask_om_cue2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 15 + col_offset])
	mask_om_p1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 17 + col_offset])
	mask_om_p2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 19 + col_offset])
	mask_prem = map(x -> x == 4,  D[:,2])

	cue[map((x,y) -> x != 0 || y != 0, 
			D[:,17+col_offset], D[:,19+col_offset])] .= 5

	RT[map((x,y) -> x == 5 && y != 0, 
			cue, D[:, 17+col_offset])] = (D[map((x,y) -> x == 5 && y != 0, 
												cue, D[:, 17+col_offset]), 18+col_offset] - 
											D[map((x,y) -> x == 5 && y != 0, 
												cue, D[:, 17+col_offset]), 17+col_offset]) / 100.0
	RT[map((x,y) -> x == 5 && y != 0, 
		cue, D[:, 19+col_offset])] = (D[map((x,y) -> x == 5 && y != 0, 
											cue, D[:, 19+col_offset]), 20+col_offset] - 
										D[map((x,y) -> x == 5 && y != 0, 
											cue, D[:, 19+col_offset]), 19+col_offset]) / 100.0

	if mod(id_number, 2) == 0 
		choice[mask_corr_cue1 .| mask_corr_p1 .| mask_incorr_cue2 .| mask_incorr_p2] .= 8
		choice[mask_corr_cue2 .| mask_corr_p2 .| mask_incorr_cue1 .| mask_incorr_p1] .= 2
		outcome[mask_corr_cue1 .| mask_corr_p1] .= 1
		outcome[mask_corr_cue2 .| mask_corr_p2] .= 4
		cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8
		cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2
		RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 16 + col_offset] - 
									D[map(x -> x == 2, cue), 15 + col_offset]) / 100.0
		RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 14 + col_offset] - 
									D[map(x -> x == 8, cue), 13 + col_offset]) / 100.0

		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 11 + col_offset])] .= 8
		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 12 + col_offset])] .= 2
	else
		choice[mask_corr_cue1 .| mask_corr_p1 .| mask_incorr_cue2 .| mask_incorr_p2] .= 2
		choice[mask_corr_cue2 .| mask_corr_p2 .| mask_incorr_cue1 .| mask_incorr_p1] .= 8
		outcome[mask_corr_cue1 .| mask_corr_p1] .= 4
		outcome[mask_corr_cue2 .| mask_corr_p2] .= 1
		cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
		cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8
		RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 14 + col_offset] - 
									D[map(x -> x == 2, cue), 13 + col_offset]) / 100.0
		RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 16 + col_offset] - 
									D[map(x -> x == 8, cue), 15 + col_offset]) / 100.0
		
		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 11 + col_offset])] .= 2
		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 12 + col_offset])] .= 8
	end

	RT[mask_om_cue1 .| mask_om_cue2 .| mask_om_p1 .| mask_om_p2] .= RT_max
	
	append!(
			df,
			DataFrame(
					ID = fill(subj_id, n_trials), 
					trial = collect(1:n_trials),
					date = fill(date, n_trials),
					session = fill("probe", n_trials),
					cue = cue, 
					choice = choice, 
					RT = RT, 
					outcome = outcome
					)
			)
end

function append_subj!(
					df, 
					D, 
					subj_id::String, 
					date, 
					::probe_1vs1_t;
					col_offset=2,
					reversed=false
					)

	RT_max = 20.0 # sec
	n_trials = size(D,1)

	outcome = zeros(Int64, n_trials)
	choice = zeros(Int64, n_trials)
	cue = zeros(Int64, n_trials)
	RT = zeros(Float64, n_trials)
	
	id_number_idx = findlast(isequal('_'), subj_id)
	id_number = tryparse(Int64, subj_id[id_number_idx+1:end])

	# p1 : probe cue playing and route 1 lever (left) was set as correct
	# p2 : probe cue playing and route 2 lever (right) was set as correct

	mask_corr_cue1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,13 + col_offset])
	mask_corr_cue2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,15 + col_offset])
	mask_incorr_cue1 = map((x,y) -> x == 1 && y != 0, D[:,2], D[:,13 + col_offset])
	mask_incorr_cue2 = map((x,y) -> x == 3 && y != 0, D[:,2], D[:,15 + col_offset])
	mask_corr_p1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,17 + col_offset])
	mask_corr_p2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,19 + col_offset])
	mask_incorr_p1 = map((x,y) -> x == 1 && y != 0, D[:,2], D[:,17 + col_offset])
	mask_incorr_p2 = map((x,y) -> x == 3 && y != 0, D[:,2], D[:,19 + col_offset])
	mask_om_cue1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 13 + col_offset])
	mask_om_cue2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 15 + col_offset])
	mask_om_p1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 17 + col_offset])
	mask_om_p2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 19 + col_offset])
	mask_prem = map(x -> x == 4,  D[:,2])

	cue[map((x,y) -> x != 0 || y != 0, 
			D[:,17+col_offset], D[:,19+col_offset])] .= 5

	RT[map((x,y) -> x == 5 && y != 0, 
			cue, D[:, 17+col_offset])] = (D[map((x,y) -> x == 5 && y != 0, 
												cue, D[:, 17+col_offset]), 18+col_offset] - 
											D[map((x,y) -> x == 5 && y != 0, 
												cue, D[:, 17+col_offset]), 17+col_offset]) / 100.0
	RT[map((x,y) -> x == 5 && y != 0, 
		cue, D[:, 19+col_offset])] = (D[map((x,y) -> x == 5 && y != 0, 
											cue, D[:, 19+col_offset]), 20+col_offset] - 
										D[map((x,y) -> x == 5 && y != 0, 
											cue, D[:, 19+col_offset]), 19+col_offset]) / 100.0

	outcome[mask_corr_cue1 .| mask_corr_p1] .= 1
	outcome[mask_corr_cue2 .| mask_corr_p2] .= 1

	if mod(id_number, 2) == 0 
		choice[mask_corr_cue1 .| mask_corr_p1 .| mask_incorr_cue2 .| mask_incorr_p2] .= 8
		choice[mask_corr_cue2 .| mask_corr_p2 .| mask_incorr_cue1 .| mask_incorr_p1] .= 2
		
		cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8
		cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2
		RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 16 + col_offset] - 
									D[map(x -> x == 2, cue), 15 + col_offset]) / 100.0
		RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 14 + col_offset] - 
									D[map(x -> x == 8, cue), 13 + col_offset]) / 100.0

		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 11 + col_offset])] .= 8
		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 12 + col_offset])] .= 2
	else
		choice[mask_corr_cue1 .| mask_corr_p1 .| mask_incorr_cue2 .| mask_incorr_p2] .= 2
		choice[mask_corr_cue2 .| mask_corr_p2 .| mask_incorr_cue1 .| mask_incorr_p1] .= 8

		cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
		cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8

		RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 14 + col_offset] - 
									D[map(x -> x == 2, cue), 13 + col_offset]) / 100.0
		RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 16 + col_offset] - 
									D[map(x -> x == 8, cue), 15 + col_offset]) / 100.0
		
		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 11 + col_offset])] .= 2
		choice[map((x,y) -> x == true && y != 0, mask_prem, D[:, 12 + col_offset])] .= 8
	end

	RT[mask_om_cue1 .| mask_om_cue2 .| mask_om_p1 .| mask_om_p2] .= RT_max
	
	append!(
			df,
			DataFrame(
					ID = fill(subj_id, n_trials), 
					trial = collect(1:n_trials),
					date = fill(date, n_trials),
					session = fill("probe", n_trials),
					cue = cue, 
					choice = choice, 
					RT = RT, 
					outcome = outcome
					)
			)
end

function append_subj!(
					df, 
					D, 
					subj_id::String, 
					date, 
					::probe_compound_1vs1_t;
					col_offset=0,
					reversed=false
					)
	
	RT_max = 8.0 # sec
	n_trials = size(D,1)

	outcome = zeros(Int64, size(D,1))
	choice = zeros(Int64, size(D,1))
	cue = zeros(Int64, size(D,1))
	RT = zeros(Float64, size(D,1))

	id_number_idx = findlast(isequal('_'), subj_id)
	id_number = tryparse(Int64, subj_id[id_number_idx + 1 : end])

	mask_corr_cue1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,17])
	mask_corr_cue2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,21])
	mask_incorr_cue1 = map((x,y) -> (x == 1 || x == 3) && y != 0, D[:,2], D[:,17])
	mask_incorr_cue2 = map((x,y) -> (x == 1 || x == 3) && y != 0, D[:,2], D[:,21])
	mask_om_cue1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 17])
	mask_om_cue2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 21])
	mask_prem = map(x -> x == 4,  D[:,2])

	mask_choice_p1 = map((x,y) -> x != 0 && y != 0, D[:, 23], D[:, 29])
	mask_choice_p2 = map((x,y) -> x != 0 && y != 0, D[:, 23], D[:, 31])
	mask_om_p = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 23])

	outcome[mask_corr_cue1] .= 1
	outcome[mask_corr_cue2] .= 1
	outcome[map((x,y) -> x != 0 || y != 0, D[:, 30], D[:, 32])] .= 1

	if reversed
		if mod(id_number, 2) != 0 # traditionally mod(id_number, 2) == 0
			choice[mask_corr_cue1 .| mask_incorr_cue2 .| mask_choice_p1] .= 8 
			choice[mask_corr_cue2 .| mask_incorr_cue1 .| mask_choice_p2] .= 2 

			choice[map(x -> x != 0, D[:, 13])] .= 8
			choice[map(x -> x != 0, D[:, 14])] .= 2

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8 
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2
			cue[mask_choice_p1 .| mask_choice_p2 .| mask_om_p] .= 5 

			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 22] - 
										D[map(x -> x == 2, cue), 21]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 18] - 
										D[map(x -> x == 8, cue), 17]) / 100.0
			RT[map(x -> x == 5, cue)] = (	D[map(x -> x == 5, cue), 24] - 
										D[map(x -> x == 5, cue), 23]) / 100.0

		else
			choice[mask_corr_cue1 .| mask_incorr_cue2 .| mask_choice_p1] .= 2 
			choice[mask_corr_cue2 .| mask_incorr_cue1 .| mask_choice_p2] .= 8 

			choice[map(x -> x != 0, D[:, 13])] .= 2
			choice[map(x -> x != 0, D[:, 14])] .= 8

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8
			cue[mask_choice_p1 .| mask_choice_p2 .| mask_om_p] .= 5 

			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 18] - 
										D[map(x -> x == 2, cue), 17]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 22] - 
										D[map(x -> x == 8, cue), 21]) / 100.0
			RT[map(x -> x == 5, cue)] = (D[map(x -> x == 5, cue), 24] - 
										D[map(x -> x == 5, cue), 23]) / 100.0
		end
	else
		if mod(id_number, 2) == 0
			choice[mask_corr_cue1 .| mask_incorr_cue2 .| mask_choice_p1] .= 8 
			choice[mask_corr_cue2 .| mask_incorr_cue1 .| mask_choice_p2] .= 2 

			choice[map(x -> x != 0, D[:, 13])] .= 8
			choice[map(x -> x != 0, D[:, 14])] .= 2

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8 
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2
			cue[mask_choice_p1 .| mask_choice_p2 .| mask_om_p] .= 5 

			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 22] - 
										D[map(x -> x == 2, cue), 21]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 18] - 
										D[map(x -> x == 8, cue), 17]) / 100.0
			RT[map(x -> x == 5, cue)] = (D[map(x -> x == 5, cue), 24] - 
										D[map(x -> x == 5, cue), 23]) / 100.0

		else
			choice[mask_corr_cue1 .| mask_incorr_cue2 .| mask_choice_p1] .= 2 
			choice[mask_corr_cue2 .| mask_incorr_cue1 .| mask_choice_p2] .= 8 

			choice[map(x -> x != 0, D[:, 13])] .= 2
			choice[map(x -> x != 0, D[:, 14])] .= 8

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8
			cue[mask_choice_p1 .| mask_choice_p2 .| mask_om_p] .= 5 

			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 18] - 
										D[map(x -> x == 2, cue), 17]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 22] - 
										D[map(x -> x == 8, cue), 21]) / 100.0
			RT[map(x -> x == 5, cue)] = (D[map(x -> x == 5, cue), 24] - 
										D[map(x -> x == 5, cue), 23]) / 100.0
		end
	end

	RT[mask_om_cue1 .| mask_om_cue2 .| mask_om_p] .= RT_max
	
	append!(
			df,
			DataFrame(
					ID = fill(subj_id, n_trials), 
					trial = collect(1:n_trials),
					date = fill(date, n_trials),
					session = fill("probe", n_trials),
					cue = cue, 
					choice = choice, 
					RT = RT, 
					outcome = outcome
					)
			)

end

function append_subj!(
					df, 
					D, 
					subj_id::String, 
					date, 
					t::Union{train_1vs1_t, train_4vs1_t}; 
					col_offset=0,
					reversed=false
					)

	RT_max = 20.0 # sec
	n_trials = size(D,1)

	outcome = zeros(Int64, n_trials)
	choice = zeros(Int64, n_trials)
	cue = zeros(Int64, n_trials)
	RT = zeros(Float64, n_trials)

	id_number_idx = findlast(isequal('_'), subj_id)
	id_number = tryparse(Int64, subj_id[id_number_idx + 1 : end])

	mask_corr_cue1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,15 + col_offset])
	mask_corr_cue2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,17 + col_offset])
	mask_incorr_cue1 = map((x,y) -> x == 1 && y != 0, D[:,2], D[:,15 + col_offset])
	mask_incorr_cue2 = map((x,y) -> x == 3 && y != 0, D[:,2], D[:,17 + col_offset])
	mask_om_cue1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 15 + col_offset])
	mask_om_cue2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 17 + col_offset])
	mask_prem = map(x -> x == 4,  D[:,2])

	if mod(id_number, 2) == 0 
		choice[mask_corr_cue1 .| mask_incorr_cue2] .= 8
		choice[mask_corr_cue2 .| mask_incorr_cue1] .= 2

		outcome[mask_corr_cue1] .= 1
		if typeof(t) == train_1vs1_t			
			outcome[mask_corr_cue2] .= 1
		else
			outcome[mask_corr_cue2] .= 4
		end

		cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8
		cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2

		RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 18] - 
									D[map(x -> x == 2, cue), 17]) / 100.0
		RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 16] - 
									D[map(x -> x == 8, cue), 15]) / 100.0
	else
		choice[mask_corr_cue1 .| mask_incorr_cue2] .= 2
		choice[mask_corr_cue2 .| mask_incorr_cue1] .= 8
		
		outcome[mask_corr_cue2] .= 1
		if typeof(t) == train_1vs1_t
			outcome[mask_corr_cue1] .= 1
		else
			outcome[mask_corr_cue1] .= 4
		end

		cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
		cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8

		RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 16] - 
									D[map(x -> x == 2, cue), 15]) / 100.0
		RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 18] - 
									D[map(x -> x == 8, cue), 17]) / 100.0
		
	end

	RT[mask_om_cue1 .| mask_om_cue2] .= RT_max

	st = typeof(t)==train_1vs1_t ? "train_1vs1" : "train_4vs1"

	append!(
			df,
			DataFrame(
					ID = fill(subj_id, n_trials), 
					trial = collect(1:n_trials),
					date = fill(date, n_trials),
					session= fill(st, n_trials),
					cue = cue, 
					choice = choice, 
					RT = RT, 
					outcome = outcome
					)
			)
end

function append_subj!(
					df, 
					D, 
					subj_id::String, 
					date, 
					t::train_1vs1_compound_t; 
					col_offset=0,
					reversed=false
					)
	
	RT_max = 8.0 # sec
	n_trials = size(D,1)

	outcome = zeros(Int64, size(D,1))
	choice = zeros(Int64, size(D,1))
	cue = zeros(Int64, size(D,1))
	RT = zeros(Float64, size(D,1))

	id_number_idx = findlast(isequal('_'), subj_id)
	id_number = tryparse(Int64, subj_id[id_number_idx + 1 : end])

	mask_corr_cue1 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,23])
	mask_corr_cue2 = map((x,y) -> x == 0 && y != 0, D[:,2], D[:,29])
	mask_incorr_cue1 = map((x,y) -> (x == 1 || x == 3) && y != 0, D[:,2], D[:,23])
	mask_incorr_cue2 = map((x,y) -> (x == 1 || x == 3) && y != 0, D[:,2], D[:,29])
	mask_om_cue1 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 23])
	mask_om_cue2 = map((x,y) -> x == 2 && y != 0, D[:,2], D[:, 29])
	mask_prem = map(x -> x == 4,  D[:,2])

	outcome[mask_corr_cue1] .= 1
	outcome[mask_corr_cue2] .= 1

	if reversed
		if mod(id_number, 2) != 0

			choice[mask_corr_cue1 .| mask_incorr_cue2] .= 8 
			choice[mask_corr_cue2 .| mask_incorr_cue1] .= 2 

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8 
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2
			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 30] - 
										D[map(x -> x == 2, cue), 29]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 24] - 
										D[map(x -> x == 8, cue), 23]) / 100.0

		else
			choice[mask_corr_cue1 .| mask_incorr_cue2] .= 2
			choice[mask_corr_cue2 .| mask_incorr_cue1] .= 8

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8
			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 24] - 
										D[map(x -> x == 2, cue), 23]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 30] - 
										D[map(x -> x == 8, cue), 29]) / 100.0
			
		end
	else
		if mod(id_number, 2) == 0

			choice[mask_corr_cue1 .| mask_incorr_cue2] .= 8 
			choice[mask_corr_cue2 .| mask_incorr_cue1] .= 2 

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 8 
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 2
			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 30] - 
										D[map(x -> x == 2, cue), 29]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 24] - 
										D[map(x -> x == 8, cue), 23]) / 100.0

		else
			choice[mask_corr_cue1 .| mask_incorr_cue2] .= 2
			choice[mask_corr_cue2 .| mask_incorr_cue1] .= 8

			cue[mask_corr_cue1 .| mask_incorr_cue1 .| mask_om_cue1] .= 2
			cue[mask_corr_cue2 .| mask_incorr_cue2 .| mask_om_cue2] .= 8
			RT[map(x -> x == 2, cue)] = (D[map(x -> x == 2, cue), 24] - 
										D[map(x -> x == 2, cue), 23]) / 100.0
			RT[map(x -> x == 8, cue)] = (D[map(x -> x == 8, cue), 30] - 
										D[map(x -> x == 8, cue), 29]) / 100.0
			
		end
	end

	RT[mask_om_cue1 .| mask_om_cue2] .= RT_max

	append!(
			df,
			DataFrame(
					ID = fill(subj_id, n_trials), 
					trial = collect(1:n_trials),
					date = fill(date, n_trials),
					session= fill("train_1vs1", n_trials),
					cue = cue, 
					choice = choice, 
					RT = RT, 
					outcome = outcome
					)
			)
end
