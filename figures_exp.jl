
function figure_drug(batch_ID, study, sessions, session_labels, cue_labels)

	fig_sz_inch = (6.4, 8)
	font_sz = 11

	dir = string("./exp/", batch_ID, "/", study,"/")

	df = read_JBT(dir)

	df_cb = CSV.File(string(dir,"counterbalance/", study, "_counterbalance.csv")) |> DataFrame

	filter!(row -> row.session == "probe", df)

	df.session = map(row -> df_cb[df_cb.ID.==row.ID, row.date][1], eachrow(df))

	cues = [2, 5, 8]

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[2, 2] = GridLayout()
	gd = f[3, 1] = GridLayout()
	ge = f[3, 2] = GridLayout()

	ax = [
		Axis(ga[1,1],
			xticks=(1:length(sessions[sessions.!="vehicle"]), session_labels[sessions.!="vehicle"]),
			yticks=[-1.0,-0.5,0.0,0.5,1.0],
			ylabel="Difference in CBI"),
		Axis(gb[1,1],
			yticks=0.0:25.0:100.0,
			xlabel="Cue",
			xticks=(1:length(cues), cue_labels),
			ylabel="HA responses [%]"),
		Axis(gc[1,1],
			xlabel="Cue",
			xticks=(1:length(cues), cue_labels),
			ylabel="Response time [sec]"),
		Axis(gd[1,1],
			xlabel="Cue",
			xticks=(1:length(cues), cue_labels),
			ylabel="Omissions [%]"),
		Axis(ge[1,1],
			xticks=(1:length(sessions), session_labels),
			ylabel="Premature responses [%]")
		]

	plot_diff_CBI!(ax[1], df, sessions, session_labels)
	plot_responses!(ax[2], df, 2, cues, sessions, session_labels)
	plot_RT!(ax[3], df, cues, sessions, session_labels)
	plot_omissions!(ax[4], df, cues, sessions, session_labels)
	plot_prematures!(ax[5], df, sessions, session_labels)

	xlims!(ax[1], [0.5, length(sessions[sessions.!="vehicle"])+0.5])
	xlims!(ax[5], 0.5, length(sessions)+0.5)
	ylims!(ax[1], -1.6, 1.6)

	colgap!(f.layout, Relative(0.04))
	rowgap!(f.layout, Relative(0.005))

	f[1, 2] = Legend(f, ax[2], framevisible=false, unique=true, tellwidth=false, tellheight=false)

	for (label, layout) in zip(["A", "B", "C", "D", "E"], [ga, gb, gc, gd, ge])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save(string("./figures/exp/", batch_ID, "_", study, ".eps"), f, pt_per_unit = 1)
end

function figure_probe(batch_ID, study, sessions, session_labels, cue_labels;
					ID_excluded=[], S_excluded=[], reversed=false)

	fig_sz_inch = (6.4, 8)
	font_sz = 11

	dir = string("./exp/", batch_ID, "/", study,"/")

	df = read_JBT(dir; reversed=reversed)

	filter!(row -> row.session == "probe", df)

	df_cb = CSV.File(string(dir,"counterbalance/", study, "_counterbalance.csv")) |> DataFrame
	df.session = map(row -> df_cb[df_cb.ID.==row.ID, row.date][1], eachrow(df))

	filter!(row -> !in(row.ID, ID_excluded) && !in(row.session, S_excluded), df)

	cues = [2, 5, 8]

	f = Figure(resolution = 72 .* fig_sz_inch, fontsize=font_sz)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[2, 2] = GridLayout()
	gd = f[3, 1] = GridLayout()
	ge = f[3, 2] = GridLayout()

	ax = [
		Axis(ga[1,1],
			xticks=(1:length(sessions), session_labels),
			ylabel="CBI"),
		Axis(gb[1,1],
			yticks=0.0:25.0:100.0,
			xlabel="Cue",
			xticks=(1:length(cues), cue_labels),
			ylabel="HA responses [%]"),
		Axis(gc[1,1],
			xlabel="Cue",
			xticks=(1:length(cues), cue_labels),
			ylabel="Response time [sec]"),
		Axis(gd[1,1],
			xlabel="Cue",
			xticks=(1:length(cues), cue_labels),
			ylabel="Omissions [%]"),
		Axis(ge[1,1],
			xticks=(1:length(sessions), session_labels),
			ylabel="Premature responses [%]")
		]

	plot_CBI!(ax[1], df, sessions, session_labels)
	plot_responses!(ax[2], df, 2, cues, sessions, session_labels)
	plot_RT!(ax[3], df, cues, sessions, session_labels)
	plot_omissions!(ax[4], df, cues, sessions, session_labels)
	plot_prematures!(ax[5], df, sessions, session_labels)

	xlims!(ax[1], 0.5, length(sessions)+0.5)
	xlims!(ax[5], 0.5, length(sessions)+0.5)
	ylims!(ax[1], -1.4, 1.4)

	colgap!(f.layout, Relative(0.04))
	rowgap!(f.layout, Relative(0.005))

	f[1, 2] = Legend(f, ax[2], framevisible=false, unique=true, tellwidth=false, tellheight=false)

	for (label, layout) in zip(["A", "B", "C", "D", "E"], [ga, gb, gc, gd, ge])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
	        halign = :right)
	end

	save(string("./figures/exp/", batch_ID, "_", study, ".eps"), f, pt_per_unit = 1)
end

figure_probe("HO2",
			"1vs1",
			["1vs1_2"],
			["Session 2"],
			["2 kHz", "5 kHz", "8 kHz"];
			ID_excluded=["HO2_4", "HO2_11"],
			S_excluded=["1vs1_1"]
			)

figure_probe("HO2",
			"4vs1",
			["4vs1_1", "4vs1_2"],
			["Session 1", "Session 2"],
			["HC", "AC", "LC"]
			)

figure_drug("HO2",
			"24hket",
			["vehicle", "ketamine_24h"],
			["Vehicle", "Ketamine (24h)"],
			["HC", "AC", "LC"]
			)

figure_drug("HO2",
			"amph",
			["vehicle", "amphetamine"],
			["Vehicle", "Amphetamine"],
			["HC", "AC", "LC"]
			)

figure_drug("HO2",
			"ket",
			["vehicle", "ketamine"],
			["Vehicle", "Ketamine"],
			["HC", "AC", "LC"]
			)

figure_drug("HO2",
			"ket_ds",
			["vehicle", "ketamine_03", "ketamine_1", "ketamine_3"],
			["V", "K 0.3", "K 1.0", "K 3.0"],
			["HC", "AC", "LC"]
			)

figure_probe("HO3",
			"1vs1",
			["1vs1_1", "1vs1_2"],
			["Session 1", "Session 2"],
			["2 kHz\n+ LL", "LL", "8 kHz\n+ LL"]
			)

figure_drug("HO3",
			"ket_1vs1",
			["vehicle", "ketamine"],
			["Vehicle", "Ketamine"],
			["2 kHz\n+ LL", "LL", "8 kHz\n+ LL"]
			)

figure_probe("HO3",
			"1vs1_reversed",
			["1vs1_reversed_1", "1vs1_reversed_2"],
			["Session 1", "Session 2"],
			["2 kHz\n+ LL", "LL", "8 kHz\n+ LL"];
			reversed=true
			)
