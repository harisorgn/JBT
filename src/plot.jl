
function plot_responses!(ax, df, choice, cues, sessions, session_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	k = length(sessions) > 1 ? LinRange(-0.08*length(sessions), 0.08*length(sessions), length(sessions)) : [0.0]

	d = Uniform(-0.05, 0.05)

	ID = unique(df.ID)

	for (s,j) in zip(sessions, eachindex(sessions))

		Y = [responses(df, i, s, choice, c) for i in ID, c in cues]
		x = (1:length(cues)) .+ k[j]

		errorbars!(ax, x, nanmean(Y,1)[:], nanstd(Y,1)[:]/sqrt(length(ID)),
					color=colormap[j], whiskerwidth=10, label=session_labels[j])

		scatter!(ax, x, nanmean(Y,1)[:],
				color=colormap[j], marker=:diamond, markersize=12)

		scatter!(ax, repeat(x, inner=length(ID)) + rand(d, length(ID)*length(x)), Y[:],
				color=(colormap[j], 0.2), markersize=6)

	end
end

function plot_RT!(ax, df, cues, sessions, session_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	k = length(sessions) > 1 ? LinRange(-0.08*length(sessions), 0.08*length(sessions), length(sessions)) : [0.0]

	d = Uniform(-0.05, 0.05)

	ID = unique(df.ID)

	for (s,j) in zip(sessions, eachindex(sessions))

		Y = [RT(df, i, s, c, c) for i in ID, c in cues]
		x = (1:length(cues)) .+ k[j]

		errorbars!(ax, x, nanmean(Y,1)[:], nanstd(Y,1)[:]/sqrt(length(ID)),
					color=colormap[j], whiskerwidth=10, label=session_labels[j])

		scatter!(ax, x, nanmean(Y,1)[:],
				color=colormap[j], marker=:diamond, markersize=12)

		scatter!(ax, repeat(x, inner=length(ID)) + rand(d, length(ID)*length(x)), Y[:],
				color=(colormap[j], 0.2), markersize=6)

	end

end

function plot_omissions!(ax, df, cues, sessions, session_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	k = length(sessions) > 1 ? LinRange(-0.08*length(sessions), 0.08*length(sessions), length(sessions)) : [0.0]

	d = Uniform(-0.05, 0.05)

	ID = unique(df.ID)

	for (s,j) in zip(sessions, eachindex(sessions))

		Y = [omissions(df, i, s, c) for i in ID, c in cues]
		x = (1:length(cues)) .+ k[j]

		errorbars!(ax, x, nanmean(Y,1)[:], nanstd(Y,1)[:]/sqrt(length(ID)),
					color=colormap[j], whiskerwidth=10, label=session_labels[j])

		scatter!(ax, x, nanmean(Y,1)[:],
				color=colormap[j], marker=:diamond, markersize=12)

		scatter!(ax, repeat(x, inner=length(ID)) + rand(d, length(ID)*length(x)), Y[:],
				color=(colormap[j], 0.2), markersize=6)

	end
end

function plot_prematures!(ax, df, sessions, session_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Uniform(-0.1, 0.1)

	ID = unique(df.ID)

	Y = [prematures(df, i, s) for i in ID, s in sessions]

	for i=1:length(sessions)

		errorbars!(ax, [i], [nanmean(Y,1)[i]], [nanstd(Y,1)[i]/sqrt(length(ID))],
					whiskerwidth=10, color=colormap[i])

		scatter!(ax, [i], [nanmean(Y,1)[i]],
					marker=:diamond, markersize=12, color=colormap[i])

		scatter!(ax, fill(i, length(ID)) + rand(d, length(ID)), Y[:,i],
					color=(colormap[i], 0.2), markersize=6)
	end
end

function plot_diff_CBI!(ax, df, sessions, session_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Uniform(-0.1, 0.1)

	ID = unique(df.ID)

	Y = [diff_CBI(df[(df.ID.==i) .& (df.session.==s),:],
				df[(df.ID.==i) .& (df.session.=="vehicle"),:])
		for i in ID, s in sessions[sessions.!="vehicle"]]

	for i=1:length(sessions[sessions.!="vehicle"])

		tst = OneSampleTTest(Y[:,i])
		pv = pvalue(tst)

		if pv < 0.001
			scatter!(ax, [i-0.1,i,i+0.1], fill(1.35,3),
					marker='*', markersize=15, color=:black)
		elseif pv < 0.01
			scatter!(ax, [i-0.05,i+0.05], fill(1.35,2),
					marker='*', markersize=15, color=:black)
		elseif pv < 0.05
			scatter!(ax, [i], fill(1.35,1),
					marker='*', markersize=15, color=:black)
		end

		errorbars!(ax, [i], [nanmean(Y,1)[i]], [nanstd(Y,1)[i]/sqrt(length(ID))],
					whiskerwidth=10, color=colormap[i+1])

		scatter!(ax, [i], [nanmean(Y,1)[i]],
					marker=:diamond, markersize=12, color=colormap[i+1])

		scatter!(ax, fill(i, length(ID)) + rand(d, length(ID)), Y[:,i],
					color=(colormap[i+1], 0.2), markersize=6)
	end
end

function plot_CBI!(ax, df, sessions, session_labels)

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Uniform(-0.1, 0.1)

	ID = unique(df.ID)

	Y = [CBI(df[(df.ID.==i) .& (df.session.==s),:]) for i in ID, s in sessions]

	for i=1:length(sessions)

		tst = OneSampleTTest(Y[:,i])
		pv = pvalue(tst)

		if pv < 0.001
			scatter!(ax, [i-0.1,i,i+0.1], fill(1.2,3),
					marker='*', markersize=15, color=:black)
		elseif pv < 0.01
			scatter!(ax, [i-0.05,i+0.05], fill(1.2,2),
					marker='*', markersize=15, color=:black)
		elseif pv < 0.05
			scatter!(ax, [i], fill(1.2,1),
					marker='*', markersize=15, color=:black)
		end

		errorbars!(ax, [i], [nanmean(Y,1)[i]], [nanstd(Y,1)[i]/sqrt(length(ID))],
					whiskerwidth=10, color=colormap[i])

		scatter!(ax, [i], [nanmean(Y,1)[i]],
					marker=:diamond, markersize=12, color=colormap[i])

		scatter!(ax, fill(i, length(ID)) + rand(d, length(ID)), Y[:,i],
					color=(colormap[i], 0.2), markersize=6)
	end
end
