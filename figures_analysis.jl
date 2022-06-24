using JBT
using CairoMakie
using ColorSchemes
using Distributions: Uniform
using Serialization
using ArviZ: hdi
using KernelDensity: kde

function win_stay_lose_shift(df)

	fig_sz_inch = (6.4, 4.8)
	font_sz = 12
	colormap = ColorSchemes.seaborn_colorblind.colors

	df_wsls = DataFrame(
					choice = Int64[],
					H4 = Float64[],
					H0 = Float64[],
					L1 = Float64[],
					L0 = Float64[]
					)

	for df_date in groupby(df, :date)

		d = combine(groupby(df_date,:ID), AsTable([:choice, :cue, :outcome]) =>
										(t -> [(
												wsls(t,2,5,2,4),
												wsls(t,2,5,2,0),
												wsls(t,2,5,8,1),
												wsls(t,2,5,8,0)
											)]) =>
										[:H4, :H0, :L1, :L0])

		d.choice = fill(2, length(d.ID))
		append!(df_wsls, d[!, Not(:ID)])

		d = combine(groupby(df_date,:ID), AsTable([:choice, :cue, :outcome]) =>
										(t -> [(
												wsls(t,8,5,2,4),
												wsls(t,8,5,2,0),
												wsls(t,8,5,8,1),
												wsls(t,8,5,8,0)
											)]) =>
										[:H4, :H0, :L1, :L0])

		d.choice = fill(8, length(d.ID))
		append!(df_wsls, d[!, Not(:ID)])
	end

	filter!(isfinite ∘ sum, df_wsls)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)
	ax = Axis(f[1,1],
			xticks=(1:4, ["A₋₁=HA,\nR₋₁=HR", "A₋₁=HA,\nR₋₁=0", "A₋₁=LA,\nR₋₁=LR", "A₋₁=LA,\nR₋₁=0"]),
			ylabel="Conditional probability ratio"
			)

	d = df_wsls[df_wsls.choice.==2, Not(:choice)]
	for (c,i) in zip(eachcol(d), 1:ncol(d))
		violin!(ax, fill(i, length(c)), c,
				datalimits=extrema, color=colormap[1], side=:right, label="High-reward action (HA)")
	end

	d = df_wsls[df_wsls.choice.==8, Not(:choice)]
	for (c,i) in zip(eachcol(d), 1:ncol(d))
		violin!(ax, fill(i, length(c)), c,
				datalimits=extrema, color=colormap[2], side=:left, label="Low-reward action (LA)")
	end

	hlines!(ax,[1.0], color=:black, linestyle=:dash, linewidth=0.5)
	ylims!(ax, -0.5, 6)

	axislegend("Current action:", unique = true, patchsize = (10, 10), labelsize=9, titlesize=9, position=:lt)

	save("./figures/analysis/wsls.eps",f, pt_per_unit=1)
end

function CBI_over_trials(df; n_blocks=4)

	fig_sz_inch = (6.4, 4)
	font_sz = 11
	colormap = ColorSchemes.seaborn_colorblind.colors

	df_CBI = DataFrame(["B$(i)" => Float64[] for i=1:n_blocks])

	for df_date in groupby(df, :date)
		d = combine(groupby(df_date,:ID), AsTable([:choice, :cue, :trial]) =>
										(x -> [CBI_blocks(x, n_blocks=n_blocks)]) =>
										["B$(i)" for i=1:n_blocks])
		append!(df_CBI, d[!, Not(:ID)])
	end

	filter!(isfinite ∘ sum, df_CBI)

	d = Uniform(-0.1, 0.1)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[1, 2] = GridLayout()

	ax = [
		Axis(ga[1,1],
			xticks=(1:5, ["B$(i)" for i=1:n_blocks]),
			ylabel="CBI"),
		Axis(gb[1,1],
			xticks=(1:4, ["B1→B2", "B2→B3", "B3→B4", "B4→B5"]),
			ylabel="Difference in CBI")
		]

	for (c,i) in zip(eachcol(df_CBI), 1:ncol(df_CBI))

		scatter!(ax[1], rand(d,length(c)).+i, c, marker=:circle, markersize=6, color=(colormap[1],0.05))

		errorbars!(ax[1], [i], [nanmean(c)], [nanstd(c)/sqrt(length(c))],
					whiskerwidth=8, color=:black)

		scatter!(ax[1], [i], [nanmean(c)],
					marker=:diamond, markersize=8, color=:black)
	end

	hlines!(ax[1],[0.0], color=:black, linestyle=:dash, linewidth=0.5)

	for (c, c_prev, i) in zip(eachcol(df_CBI)[2:end], eachcol(df_CBI)[1:end-1], 1:ncol(df_CBI)-1)

		scatter!(ax[2], rand(d,length(c)).+i, c - c_prev, marker=:circle, markersize=6, color=(colormap[1],0.1))

		errorbars!(ax[2], [i], [nanmean(c) - nanmean(c_prev)], [(nanstd(c)+nanstd(c_prev))/sqrt(length(c))],
					whiskerwidth=8, color=:black)

		scatter!(ax[2], [i], [nanmean(c) - nanmean(c_prev)],
					marker=:diamond, markersize=8, color=:black)
	end

	hlines!(ax[2],[0.0], color=:black, linestyle=:dash, linewidth=0.5)

	colsize!(f.layout, 1, Auto(1.2))
	colsize!(f.layout, 2, Auto(1.2))
	colgap!(f.layout, Relative(0.04))

	for (label, layout) in zip(["A", "B"], [ga, gb])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save("./figures/analysis/CBI_blocks.eps",f, pt_per_unit=1)
end

function group_posterior_amph()

	fig_sz_inch = (6.4,9)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind.colors

	df = get_cohort_df(["amph"])
	(choices, data) = df_to_JBT(df, cue_reward_rate)
	model = logistic_past(choices, data)
	chain = deserialize("chain_amph_learning.jls")

	chain_prior = sample(model, Prior(), 4000)
	df_prior = DataFrame(chain_prior)
	df_c = DataFrame(chain)

	n_samples = nrow(df_c)

	d_condition = condition_dict(df)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[1, 2] = GridLayout()
	gc = f[2, 1] = GridLayout()
	gd = f[2, 2] = GridLayout()
	ge = f[3, 1] = GridLayout()
	gf = f[3, 2] = GridLayout()

	ax = [
		Axis(ga[1,1],
			xticks=(0:2, ["Prior", "Vehicle", "Amph"]),
			ylabel=L"Population bias, $\mu$"
			),
		Axis(gb[1,1],
			xlabel=L"Difference in $\mu$"
			),
		Axis(gc[1,1],
			xticks=(0:2, ["Prior", "Vehicle", "Amph"]),
			ylabel=L"Population effect of previous trial, $\xi$"
			),
		Axis(gd[1,1],
			xlabel=L"Difference in $\xi$"
			),
		Axis(ge[1,1],
			xticks=(0:2, ["Prior", "Vehicle", "Amph"]),
			ylabel=L"Population effect of learning, $\eta$"
			),
		Axis(gf[1,1],
			xlabel=L"Difference in $\eta$"
			)
		]

	condition_v = ["vehicle", "amphetamine"]
	param_v = ["μ_b", "μ_p", "μ_r"]

	for (param,i) in zip(param_v, 1:2:6)
		c = group(chain, param)

		violin!(ax[i], fill(0, length(df_prior[!, "μ_b[1]"])), df_prior[!, "μ_b[1]"],
				side=:right, color=(:slategray, 0.4))

		for (condition,j) in zip(condition_v, eachindex(condition_v))
			violin!(ax[i], fill(j, n_samples), c[:,d_condition[condition],:][:].data,
				side=:right, color=colormap[j])
		end

		density!(ax[i+1],
				c[:, d_condition["amphetamine"],:][:].data - c[:,d_condition["vehicle"],:][:].data,
				color=colormap[2])

		vlines!(ax[i+1],[0.0], color=:black, linestyle=:dash, linewidth=0.8)
		xlims!(ax[i+1], -4.0, 4.0)
	end

	hidexdecorations!.(ax[1:4], grid=false, label=false)
	colgap!(f.layout, Relative(0.02))
	rowgap!(f.layout, Relative(0.0))
	colsize!(f.layout, 1, Auto(1.2))
	rowsize!(f.layout, 1, Auto(2))
	rowsize!(f.layout, 2, Auto(2))
	rowsize!(f.layout, 3, Auto(2))

	for (label, layout) in zip(["A", "B", "C", "D", "E", "F"], [ga, gb, gc, gd, ge, gf])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save("./figures/analysis/group_post_amph.eps", f, pt_per_unit = 1)
end

function subject_posterior_lapse_amph()

	fig_sz_inch = (6.4, 4.8)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind6.colors

	df = get_cohort_df(["amph"])

	(choices, data) = df_to_JBT(df, cue_reward_rate)
	model = logistic_past(choices, data)
	chain = deserialize("chain_amph_learning.jls")

	chain_prior = sample(model, Prior(), 4000)
	df_prior = DataFrame(chain_prior)
	df_c = DataFrame(chain)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)
	ax = Axis(f[1,1],
			xticks=(0:1, ["Prior", "Posterior"]),
			ylabel=L"Subject lapse rate $\epsilon$")

	ε_prior = reduce(
					vcat,
					[df_prior[!, "μ_ε"] .+ df_prior[!, "ε_norm[$(i)]"].*df_prior[!, "σ_ε"]
					for i=1:data.n_subjects]
					)
	ε = reduce(
				vcat,
				[df_c[!, "μ_ε"] .+ df_c[!, "ε_norm[$(i)]"].*df_c[!, "σ_ε"]
				for i=1:data.n_subjects]
				)

	violin!(ax, fill(0, length(ε_prior)), cdf.(Normal(0,1), ε_prior),
			datalimits=extrema, side=:right, color=(:slategray, 0.4))
	violin!(ax, fill(1, length(ε)), cdf.(Normal(0,1), ε),
			datalimits=extrema, side=:right, color=colormap[1])

	save("./figures/analysis/subj_post_lapse_amph.eps", f, pt_per_unit = 1)
end

function subject_prior_lapse()

	fig_sz_inch = (6.4, 4.8)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind6.colors

	d = Dict(
			"./exp/probe/baseline/baseline_no_effect/"=>"baseline",
			"./exp/probe/ketamine/ketamine_no_effect/"=>"ketamine",
			"./exp/probe/ketamine/vehicle_no_effect/"=>"vehicle"
			)
	df = get_batch_df(d)

	(choices, data) = df_to_JBT(df, cue_reward_rate)
	model = logistic_past(choices, data)

	chain_prior = sample(model, Prior(), 4000)
	df_prior = DataFrame(chain_prior)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)
	ax = Axis(f[1,1],
			xlabel=L"Subject lapse rate $\epsilon$",
			ylabel="Prior probability density")

	ε_prior = reduce(vcat, [df_prior[!, "μ_ε"] .+ df_prior[!, "ε_norm[$(i)]"].*df_prior[!, "σ_ε"] for i=1:13])

	density!(ax, cdf.(Normal(0,1), ε_prior), color=(:slategray, 0.8))
	xlims!(ax,0,1)
	save("./figures/analysis/prior_ε.eps", f, pt_per_unit = 1)
end

function group_posterior()

	fig_sz_inch = (6.4,9)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Dict(
			"./exp/probe/baseline/baseline_no_effect/"=>"baseline",
			"./exp/probe/ketamine/ketamine_no_effect/"=>"ketamine",
			"./exp/probe/ketamine/vehicle_no_effect/"=>"vehicle"
			)
	df = get_batch_df(d)
	(choices, data) = df_to_JBT(df, cue_reward_rate)
	model = logistic_past(choices, data)
	chain_prior = sample(model, Prior(), 4000)
	df_prior = DataFrame(chain_prior)

	chain_e = deserialize("chain_learning_effect_batch.jls")
	df_e = DataFrame(chain_e)

	chain_ne = deserialize("chain_learning_no_effect_batch.jls")
	df_ne = DataFrame(chain_ne)

	d_condition = condition_dict(df)

	n_samples = nrow(df_e)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[3, 1] = GridLayout()

	ax = [
		Axis(ga[1,1],
			ylabel=L"Population bias, $\mu$",
			xticks=(0:3, ["Prior", "Baseline", "Vehicle", "Ketamine"])),
		Axis(gb[1,1],
			ylabel=L"Population effect of previous trial, $\xi$",
			xticks=(0:3, ["Prior", "Baseline", "Vehicle", "Ketamine"])),
		Axis(gc[1,1],
			ylabel=L"Population effect of learning, $\eta$",
			xticks=(0:3, ["Prior", "Baseline", "Vehicle", "Ketamine"]))
		]

	condition_v = ["baseline", "vehicle", "ketamine"]
	param_v = ["μ_b", "μ_p", "μ_r"]

	for (param,i) in zip(param_v, 1:3)

		violin!(ax[i], fill(0, length(df_prior[!, "μ_b[1]"])), df_prior[!, "μ_b[1]"],
				color=(:slategray, 0.4))

		for (condition,j) in zip(condition_v, eachindex(condition_v))
			violin!(
					ax[i],
					fill(j, n_samples),
					df_e[!, string(param, "[$(d_condition[condition])]")],
					side=:right,
					color=colormap[1],
					label="Original studies"
					)
			violin!(
					ax[i],
					fill(j, n_samples),
					df_ne[!, string(param, "[$(d_condition[condition])]")],
					side=:left,
					color=colormap[2],
					label="Replication studies"
					)
		end
	end

	ylims!.(ax,-2.7,2.7)
	hlines!.(ax,[0.0], color=:black, linestyle=:dash, linewidth=0.5)
	hidexdecorations!.(ax[1:2], grid=false, label=false)
	axislegend(ax[1], unique = true, labelsize=7, position=:rb, patchsize = (5,5))

	for (label, layout) in zip(["A", "B", "C"], [ga, gb, gc])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save("./figures/analysis/group_post.eps", f, pt_per_unit = 1)
end

function group_posterior_diff()

	fig_sz_inch = (6.4,9)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind.colors

	d = Dict(
			"./exp/probe/baseline/baseline_no_effect/"=>"baseline",
			"./exp/probe/ketamine/ketamine_no_effect/"=>"ketamine",
			"./exp/probe/ketamine/vehicle_no_effect/"=>"vehicle"
			)
	df = get_batch_df(d)

	chain_e = deserialize("chain_learning_effect_batch.jls")
	df_e = DataFrame(chain_e)

	chain_ne = deserialize("chain_learning_no_effect_batch.jls")
	df_ne = DataFrame(chain_ne)

	d_condition = condition_dict(df)

	n_samples = nrow(df_e)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[3, 1] = GridLayout()

	ax = [
		Axis(ga[1,1],
			ylabel=L"Difference in $\mu$",
			xticks=([1], ["Posterior difference"])
			),
		Axis(gb[1,1],
			ylabel=L"Difference in $\xi$",
			xticks=([1], ["Posterior difference"])
			),
		Axis(gc[1,1],
			ylabel=L"Difference in $\eta$",
			xticks=([1], ["Posterior difference"])
			)
		]

	param_v = ["μ_b", "μ_p", "μ_r"]

	for (param,i) in zip(param_v, 1:3)

			violin!(
					ax[i],
					fill(1, n_samples),
					df_e[!, string(param, "[$(d_condition["ketamine"])]")] .-
					df_e[!, string(param, "[$(d_condition["vehicle"])]")],
					side=:right,
					color=colormap[1],
					label="Original studies"
					)

			violin!(
					ax[i],
					fill(1, n_samples),
					df_ne[!, string(param, "[$(d_condition["ketamine"])]")] .-
					df_ne[!, string(param, "[$(d_condition["vehicle"])]")],
					side=:left,
					color=colormap[2],
					label="Replication studies"
					)
	end

	#ylims!.(ax,-2.7,2.7)
	hlines!.(ax,[0.0], color=:black, linestyle=:dash, linewidth=0.5)
	hidexdecorations!.(ax[1:2], grid=false, label=false)
	axislegend(ax[1], unique = true, labelsize=7, position=:rb, patchsize = (5,5))

	for (label, layout) in zip(["A", "B", "C"], [ga, gb, gc])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save("./figures/analysis/group_post_diff.eps", f, pt_per_unit = 1)
end

function subject_posterior()

	fig_sz_inch = (6.4,9)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind.colors

	df = get_batch_df(
					Dict(
						"./exp/probe/baseline/baseline_effect/"=>"baseline",
						"./exp/probe/ketamine/ketamine_effect/"=>"ketamine",
						"./exp/probe/ketamine/vehicle_effect/"=>"vehicle"
						)
					)
	(choices, data_e) = df_to_JBT(df, cue_reward_rate)

	df = get_batch_df(
					Dict(
						"./exp/probe/baseline/baseline_no_effect/"=>"baseline",
						"./exp/probe/ketamine/ketamine_no_effect/"=>"ketamine",
						"./exp/probe/ketamine/vehicle_no_effect/"=>"vehicle"
						)
					)
	(choices, data_ne) = df_to_JBT(df, cue_reward_rate)

	chain_e = deserialize("chain_learning_effect_batch.jls")
	df_e = DataFrame(chain_e)

	chain_ne = deserialize("chain_learning_no_effect_batch.jls")
	df_ne = DataFrame(chain_ne)

	d_condition = condition_dict(df)

	n_samples = nrow(df_e)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[3, 1] = GridLayout()

	ax = [
		Axis(ga[1,1],
			ylabel=L"Subject bias, $b$",
			xticks=(1:3, ["Baseline", "Vehicle", "Ketamine"])),
		Axis(gb[1,1],
			ylabel=L"Subject effect of previous trial, $w$",
			xticks=(1:3, ["Baseline", "Vehicle", "Ketamine"])),
		Axis(gc[1,1],
			ylabel=L"Subject effect of learning, $r$",
			xticks=(1:3, ["Baseline", "Vehicle", "Ketamine"]))
		]

	condition_v = ["baseline", "vehicle", "ketamine"]
	param_v = [
			("μ_b", "σ_b", "b_norm"),
			("μ_p", "σ_p", "w_p_norm"),
			("μ_r", "σ_r", "w_r_norm")
			]

	for (param,i) in zip(param_v, eachindex(param_v))
		for (condition,j) in zip(condition_v, eachindex(condition_v))

			E = [
				df_e[!, string(param[1], "[$(d_condition[condition])]")] .+
				df_e[!, string(param[3], "[$(d_condition[condition]),$(s)]")] .*
				df_e[!, string(param[2], "[$(d_condition[condition])]")]
				for s=1:data_e.n_subjects
				]

			NE = [
				df_ne[!, string(param[1], "[$(d_condition[condition])]")] .+
				df_ne[!, string(param[3], "[$(d_condition[condition]),$(s)]")] .*
				df_ne[!, string(param[2], "[$(d_condition[condition])]")]
				for s=1:data_ne.n_subjects
				]

			d = Uniform(j-0.3,j-0.05)

			scatter!(
					ax[i],
					rand(d, data_e.n_subjects) .+ 0.35,
					mean.(E),
					marker=:circle,
					markersize=6,
					color=(colormap[1],0.5),
					label="Original studies"
					)

			scatter!(
					ax[i],
					rand(d, data_ne.n_subjects),
					mean.(NE),
					marker=:circle,
					markersize=6 ,
					color=(colormap[2],0.5),
					label="Replication studies"
					)
		end
	end

	hlines!.(ax,[0.0], color=:black, linestyle=:dash, linewidth=0.5)
	hidexdecorations!.(ax[1:2], grid=false, label=false)

	ylims!(ax[1], -2.2, 7.0)
	axislegend(ax[1], unique = true, labelsize=7, position=:lt, patchsize = (5,5))

	for (label, layout) in zip(["A", "B", "C"], [ga, gb, gc])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save("./figures/analysis/subj_post.eps", f, pt_per_unit = 1)
end

function subject_posterior_lapse()

	fig_sz_inch = (6.4, 4.8)
	font_sz = 12

	colormap = ColorSchemes.seaborn_colorblind.colors

	df = get_batch_df(
					Dict(
						"./exp/probe/baseline/baseline_effect/"=>"baseline",
						"./exp/probe/ketamine/ketamine_effect/"=>"ketamine",
						"./exp/probe/ketamine/vehicle_effect/"=>"vehicle"
						)
					)
	(choices, data_e) = df_to_JBT(df, cue_reward_rate)

	df = get_batch_df(
					Dict(
						"./exp/probe/baseline/baseline_no_effect/"=>"baseline",
						"./exp/probe/ketamine/ketamine_no_effect/"=>"ketamine",
						"./exp/probe/ketamine/vehicle_no_effect/"=>"vehicle"
						)
					)
	(choices, data_ne) = df_to_JBT(df, cue_reward_rate)

	chain_e = deserialize("chain_learning_effect_batch.jls")
	df_e = DataFrame(chain_e)

	chain_ne = deserialize("chain_learning_no_effect_batch.jls")
	df_ne = DataFrame(chain_ne)

	model = logistic_past(choices, data_ne)
	chain_prior = sample(model, Prior(), 4000)
	df_prior = DataFrame(chain_prior)

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)
	ax = Axis(f[1,1],
			xticks=(0:1, ["Prior", "Posterior"]),
			ylabel=L"Subject lapse rate $\epsilon$")

	ε_prior = reduce(vcat, [df_prior[!, "μ_ε"] .+ df_prior[!, "ε_norm[$(i)]"].*df_prior[!, "σ_ε"] for i=1:13])
	ε_e = reduce(
				vcat,
				[df_e[!, "μ_ε"] .+ df_e[!, "ε_norm[$(s)]"].*df_e[!, "σ_ε"]
				for s=1:data_e.n_subjects]
				)
	ε_ne = reduce(
				vcat,
				[df_ne[!, "μ_ε"] .+ df_ne[!, "ε_norm[$(s)]"].*df_ne[!, "σ_ε"]
				for s=1:data_ne.n_subjects]
				)

	violin!(ax, fill(0, length(ε_prior)), cdf.(Normal(0,1), ε_prior),
			datalimits=extrema, color=(:slategray, 0.4))

	violin!(ax, fill(1, length(ε_e)), cdf.(Normal(0,1), ε_e),
			datalimits=extrema, side=:right, color=colormap[1],
			label="Original studies")

	violin!(ax, fill(1, length(ε_ne)), cdf.(Normal(0,1), ε_ne),
			datalimits=extrema, side=:left, color=colormap[2],
			label="Replication studies")

	axislegend(ax, unique = true, labelsize=7, position=:lt, patchsize = (5,5))

	save("./figures/analysis/subj_post_lapse.eps", f, pt_per_unit = 1)
end

df = read_JBT("./exp/probe/baseline/")

CBI_over_trials(df, n_blocks=5)
win_stay_lose_shift(df)
subject_prior_lapse()
group_posterior()
subject_posterior()
subject_posterior_lapse()
group_posterior_amph()
subject_posterior_lapse_amph()
group_posterior_diff()
