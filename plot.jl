using CairoMakie
using ColorSchemes
using LaTeXStrings
using Distributions: Uniform 
using Statistics: mean, std
using HypothesisTests

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

function RT(df, ID, session, cue)

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
					(df.choice.==cue) .& 
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

CBI(df, ID, session) = (count(x->x==2, df[(df.ID.==ID) .& (df.cue.==5) .& (df.session.==session), :choice]) -
						count(x->x==8, df[(df.ID.==ID) .& (df.cue.==5) .& (df.session.==session), :choice]))/
						(count(x->x==2, df[(df.ID.==ID) .& (df.cue.==5) .& (df.session.==session), :choice]) +
						count(x->x==8, df[(df.ID.==ID) .& (df.cue.==5) .& (df.session.==session), :choice]))

diff_CBI(df, ID, session) = CBI(df, ID, session) - CBI(df, ID, "Veh")

prematures(df, ID, session) = 100.0*count(x->x==0, df[(df.ID.==ID) .& (df.session.==session), :cue]) /
									length(df[(df.ID.==ID) .& (df.session.==session), :cue])

function plot_timeline!(ax, f, df, cue_v, S_labels)

	ID_v = unique(df.ID)
	s_v = unique(df.session)

	Y = [f(df,ID,s,t) for ID in ID_v, s in s_v, t in cue_v]

	S = reduce(vcat, [fill(i,length(ID_v)) for i in s_v])

	colormap = ColorSchemes.seaborn_colorblind.colors

	for i in eachindex(cue_v)
		
		scatter!(ax, S, Y[:,:,i][:], color=(colormap[i], 0.5), label=S_labels[i])

	end
end

function plot_responses!(ax, df, choice, C, S, S_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	k = length(S) > 1 ? LinRange(-0.08*length(S), 0.08*length(S), length(S)) : [0.0]

	d = Uniform(-0.05, 0.05)

	ID = unique(df.ID)

	for (s,j) in zip(S, eachindex(S))

		Y = [responses(df, i, s, choice, c) for i in ID, c in C]
		x = (1:length(C)) .+ k[j]

		errorbars!(ax, x, nanmean(Y,1)[:], nanstd(Y,1)[:]/sqrt(length(ID)), 
					color=colormap[j], whiskerwidth=20, label=S_labels[j])

		scatter!(ax, x, nanmean(Y,1)[:], 
				color=colormap[j], marker=:diamond, markersize=20)

		scatter!(ax, repeat(x, inner=length(ID)) + rand(d, length(ID)*length(x)), Y[:], 
				color=(colormap[j], 0.2))

	end
end

function plot_RT!(ax, df, C, S, S_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	k = length(S) > 1 ? LinRange(-0.08*length(S), 0.08*length(S), length(S)) : [0.0]
	
	d = Uniform(-0.05, 0.05)

	ID = unique(df.ID)

	for (s,j) in zip(S, eachindex(S))

		Y = [RT(df, i, s, c) for i in ID, c in C]
		x = (1:length(C)) .+ k[j]
		
		errorbars!(ax, x, nanmean(Y,1)[:], nanstd(Y,1)[:]/sqrt(length(ID)), 
					color=colormap[j], whiskerwidth=20, label=S_labels[j])

		scatter!(ax, x, nanmean(Y,1)[:], 
				color=colormap[j], marker=:diamond, markersize=20)

		scatter!(ax, repeat(x, inner=length(ID)) + rand(d, length(ID)*length(x)), Y[:], 
				color=(colormap[j], 0.2))

	end

end

function plot_omissions!(ax, df, C, S, S_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	k = length(S) > 1 ? LinRange(-0.08*length(S), 0.08*length(S), length(S)) : [0.0]
	
	d = Uniform(-0.05, 0.05)

	ID = unique(df.ID)

	for (s,j) in zip(S, eachindex(S))

		Y = [omissions(df, i, s, c) for i in ID, c in C]
		x = (1:length(C)) .+ k[j]

		errorbars!(ax, x, nanmean(Y,1)[:], nanstd(Y,1)[:]/sqrt(length(ID)), 
					color=colormap[j], whiskerwidth=20, label=S_labels[j])

		scatter!(ax, x, nanmean(Y,1)[:], 
				color=colormap[j], marker=:diamond, markersize=20)

		scatter!(ax, repeat(x, inner=length(ID)) + rand(d, length(ID)*length(x)), Y[:], 
				color=(colormap[j], 0.2))

	end
end

function plot_prematures!(ax, df, S)

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Uniform(-0.1, 0.1)

	ID = unique(df.ID)

	Y = [prematures(df, i, s) for i in ID, s in S]

	for i=1:length(S)

		errorbars!(ax, [i], [nanmean(Y,1)[i]], [nanstd(Y,1)[i]/sqrt(length(ID))], 
					whiskerwidth=20, color=colormap[i])

		scatter!(ax, [i], [nanmean(Y,1)[i]], 
					marker=:diamond, markersize=20, color=colormap[i])

		scatter!(ax, fill(i, length(ID)) + rand(d, length(ID)), Y[:,i], 
					color=(colormap[i], 0.2))
	end
end

function plot_diff_CBI!(ax, df, S)

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Uniform(-0.1, 0.1)

	ID = unique(df.ID)

	Y = [diff_CBI(df, i, s) for i in ID, s in S]

	for i=1:length(S)

		pv = pvalue(OneSampleTTest(Y[:,i]))
		@show S[i]
		@show pv

		if pv < 0.001
			scatter!(ax, [i-0.05,i,i+0.05], fill(1.35,3), 
					marker='*', markersize=25, color=:black)
		elseif pv < 0.01
			scatter!(ax, [i-0.05,i+0.05], fill(1.35,2), 
					marker='*', markersize=25, color=:black)
		elseif pv < 0.05
			scatter!(ax, [i], fill(1.35,1), 
					marker='*', markersize=25, color=:black)
		end

		errorbars!(ax, [i], [nanmean(Y,1)[i]], [nanstd(Y,1)[i]/sqrt(length(ID))], 
					whiskerwidth=20, color=colormap[i+1])

		scatter!(ax, [i], [nanmean(Y,1)[i]], 
					marker=:diamond, markersize=20, color=colormap[i+1])

		scatter!(ax, fill(i, length(ID)) + rand(d, length(ID)), Y[:,i], 
					color=(colormap[i+1], 0.2))
	end
end

function plot_CBI!(ax, df, S)

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Uniform(-0.1, 0.1)

	ID = unique(df.ID)

	Y = [CBI(df, i, s) for i in ID, s in S]

	for i=1:length(S)

		pv = pvalue(OneSampleTTest(Y[:,i]))
		@show S[i]
		@show pv

		if pv < 0.001
			scatter!(ax, [i-0.05,i,i+0.05], fill(1.2,3), 
					marker='*', markersize=25, color=:black)
		elseif pv < 0.01
			scatter!(ax, [i-0.05,i+0.05], fill(1.2,2), 
					marker='*', markersize=25, color=:black)
		elseif pv < 0.05
			scatter!(ax, [i], fill(1.2,1), 
					marker='*', markersize=25, color=:black)
		end

		errorbars!(ax, [i], [nanmean(Y,1)[i]], [nanstd(Y,1)[i]/sqrt(length(ID))], 
					whiskerwidth=20, color=colormap[i])

		scatter!(ax, [i], [nanmean(Y,1)[i]], 
					marker=:diamond, markersize=20, color=colormap[i])

		scatter!(ax, fill(i, length(ID)) + rand(d, length(ID)), Y[:,i], 
					color=(colormap[i], 0.2))
	end
end
