using Statistics: mean, std

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x;dims = y)

nanstd(x) = std(filter(!isnan,x))
nanstd(x,y) = mapslices(nanstd,x;dims = y)

responses(df, ID, session, choice, cue) = 100.0*count(
													x->x==choice, 
													df[
													(df.ID.==ID) .& 
													(df.cue.==cue) .& 
													(df.session.==session),
													:choice]
													) /
												length(
													df[
													(df.ID.==ID) .& 
													(df.choice.!=0) .&
													(df.cue.==cue) .& 
													(df.session.==session),
													:choice]
													)

function choice_prev_outcome(t, choice, cue, prev_outcome)

	ACO = map((x,y,z) -> (x==cue) && (y==choice) && (z==prev_outcome),
			t.cue[2:end],
			t.choice[2:end], 
			t.outcome[1:end-1]
			)

	CO = map((x,y) -> (x==cue) && (y==prev_outcome),
			t.cue[2:end],
			t.outcome[1:end-1]
			)

	AC = map((x,y) -> (x==cue) && (y==choice),
			t.cue[2:end],
			t.choice[2:end] 
			)

	C = count(x->x==cue, t.cue[2:end])

	#return 100.0*sum(mask)/
	#			sum((t.choice[2:end].==choice).&(t.cue[2:end].==cue))

	return (sum(ACO)/sum(CO))/(sum(AC)/C)
end

function choice_prev_choice(t, choice, cue, prev_choice)

	ACA = map((x,y,z) -> (x==cue) && (y==choice) && (z==prev_choice),
			t.cue[2:end],
			t.choice[2:end], 
			t.choice[1:end-1]
			)

	CA = map((x,y) -> (x==cue) && (y==prev_choice),
			t.cue[2:end],
			t.choice[1:end-1]
			)

	AC = map((x,y) -> (x==cue) && (y==choice),
			t.cue[2:end],
			t.choice[2:end] 
			)

	C = count(x->x==cue, t.cue[2:end])

	return (sum(ACA)/sum(CA))/(sum(AC)/C)
end

function wsls(t, choice, cue, prev_choice, prev_outcome)

	ACAO = map((x,y,z,k) -> (x==cue) && (y==choice) && (z==prev_choice) && (k==prev_outcome),
			t.cue[2:end],
			t.choice[2:end], 
			t.choice[1:end-1],
			t.outcome[1:end-1]
			)

	CAO = map((x,y,z) -> (x==cue) && (y==prev_choice) && (z==prev_outcome),
			t.cue[2:end],
			t.choice[1:end-1],
			t.outcome[1:end-1]
			)

	AC = map((x,y) -> (x==cue) && (y==choice),
			t.cue[2:end],
			t.choice[2:end] 
			)

	C = count(x->x==cue, t.cue[2:end])

	return (sum(ACAO)/sum(CAO))/(sum(AC)/C)
	#return (sum(ACAO)/sum(CAO))
end

function CBI_blocks(t; n_blocks=4)

	d = div(t.trial[end], n_blocks)
	m = mod(t.trial[end], n_blocks)

	t_blocks = [(t.trial[i*d+min(m,i)+1], t.trial[(i+1)*d+min(i+1,m)]) for i=0:(n_blocks-1)]

	b = map(t_blocks) do t_b
		CBI((choice = t.choice[(t.trial.>=t_b[1]) .& (t.trial .<=t_b[2])],
			cue = t.cue[(t.trial.>=t_b[1]) .& (t.trial .<=t_b[2])]))
	end

	return b
end

function reward_blocks(t; n_blocks=4)

	d = div(t.trial[end], n_blocks)
	m = mod(t.trial[end], n_blocks)

	t_blocks = [(t.trial[i*d+min(m,i)+1], t.trial[(i+1)*d+min(i+1,m)]) for i=0:(n_blocks-1)]

	b = map(t_blocks) do t_b
		mean(t.outcome[(t.trial.<=t_b[2])])
	end

	return b
end

function cue_reward_rate(t; window_length=12)

	trial_idx = collect(1:length(t.trial))

	R = map(trial_idx[2:end]) do trial

			mask = (t.cue.==t.cue[trial]) .& (trial_idx.<trial) .& (trial_idx.>(trial - window_length))
			all(iszero, mask) ? 0.0 : mean(t.outcome[mask])
		end
	return append!([0.0],R)
end

function reward_rate(t; window_length=12)

	trial_idx = collect(1:length(t.trial))

	R = map(trial_idx[2:end]) do trial

			mask = (trial_idx.<trial) .& (trial_idx.>(trial - window_length))
			mean(t.outcome[mask])
		end
	return append!([0.0],R)
end

function wsls_trial_index(t)

	I = map(t.choice, t.outcome) do c,o
		if c==8
			o==1 ? +1.0 : -1.0
		elseif c==2
			if any(occursin.("1vs1", t.session))
				o==1 ? -1.0 : +1.0
			else
				o==4 ? -4.0 : +4.0
			end
		end
	end

	return append!([0.0], I[1:end-1])
end

function RT(df, ID, session, choice, cue)

	if cue == 5
		return mean(
					df[
					(df.ID.==ID) .& 
					(df.cue.==cue) .& 
					(df.session.==session),
					:RT]
					)
	elseif (cue == 2) || (cue == 8)
		return mean(
					df[
					(df.ID.==ID) .& 
					(df.cue.==cue) .& 
					(df.choice.==choice) .& 
					(df.session.==session),
					:RT]
					)
	end
end

omissions(df, ID, session, cue) = 100.0*count(
											x->x==0, 
											df[
											(df.ID.==ID) .& 
											(df.cue.==cue) .& 
											(df.session.==session),
											:choice]
											) /
										length(
											df[
											(df.ID.==ID) .& 
											(df.cue.==cue) .& 
											(df.session.==session),
											:choice]
											)

CBI(t::NamedTuple) = (count(x->x==2, t.choice[t.cue.==5]) - count(x->x==8, t.choice[t.cue.==5]))/
					(count(x->x==2, t.choice[t.cue.==5]) + count(x->x==8, t.choice[t.cue.==5]))

CBI(df::DataFrame) = (count(x->x==2, df[df.cue.==5, :choice]) - count(x->x==8, df[df.cue.==5, :choice]))/
					(count(x->x==2, df[df.cue.==5, :choice]) + count(x->x==8, df[df.cue.==5, :choice]))

diff_CBI(df, df_veh) = CBI(df) - CBI(df_veh)

prematures(df, ID, session) = 100.0*count(x->x==0, df[(df.ID.==ID) .& (df.session.==session), :cue]) /
									length(df[(df.ID.==ID) .& (df.session.==session), :cue])
