
include("read.jl")
include("plot.jl")

function figure_drug(batch_ID, study, S, S_labels, C_labels)

	dir = string("./exp/", batch_ID, "/", study,"/")

	df = read_JBT(dir)

	df_cb = CSV.File(string(dir,"counterbalance/", study, "_counterbalance.csv")) |> DataFrame

	filter!(row -> row.session == "probe", df)

	df.session = map(row -> df_cb[df_cb.ID.==row.ID, row.date][1], eachrow(df))

	C = [2, 5, 8]

	f = Figure(resolution = (1280, 960), fontsize=18)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[2, 2] = GridLayout()
	gd = f[3, 1] = GridLayout()
	ge = f[3, 2] = GridLayout()

	ax = [
		Axis(ga[1,1], 
			xticks=(1:length(S[2:end]), S_labels[2:end]),
			yticks=[-1.0,-0.5,0.0,0.5,1.0], 
			ylabel="Difference in CBI"),
		Axis(gb[1,1], 
			yticks=0.0:25.0:100.0, 
			xlabel="Cue", 
			xticks=(1:length(C), C_labels), 
			ylabel="HA responses [%]"),
		Axis(gc[1,1], 
			xlabel="Cue", 
			xticks=(1:length(C), C_labels), 
			ylabel="Response time [sec]"),
		Axis(gd[1,1], 
			xlabel="Cue", 
			xticks=(1:length(C), C_labels), 
			ylabel="Omissions [%]"),
		Axis(ge[1,1], 
			xticks=(1:length(S), S_labels), 
			ylabel="Premature responses [%]")
		]

	plot_diff_CBI!(ax[1], df, S[2:end])
	plot_responses!(ax[2], df, 2, C, S, S_labels)
	plot_RT!(ax[3], df, C, S, S_labels)
	plot_omissions!(ax[4], df, C, S, S_labels)
	plot_prematures!(ax[5], df, S)

	xlims!(ax[1], [0.5, length(S[2:end])+0.5])
	ylims!(ax[1], -1.6, 1.6)

	f[1, 2] = Legend(f, ax[2], framevisible=false, unique=true, tellwidth=false, tellheight=false)

	for (label, layout) in zip(["A", "B", "C", "D", "E"], [ga, gb, gc, gd, ge])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 26,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	save(string("./figures/", batch_ID, "_", study, ".eps"),f)
end

function figure_probe(batch_ID, study, S, S_labels, C_labels; 
					ID_excluded=[], S_excluded=[], reversed=false)

	dir = string("./exp/", batch_ID, "/", study,"/")

	df = read_JBT(dir; reversed=reversed)

	filter!(row -> row.session == "probe", df)

	df_cb = CSV.File(string(dir,"counterbalance/", study, "_counterbalance.csv")) |> DataFrame
	df.session = map(row -> df_cb[df_cb.ID.==row.ID, row.date][1], eachrow(df))

	filter!(row -> !in(row.ID, ID_excluded) && !in(row.session, S_excluded), df)

	C = [2, 5, 8]

	f = Figure(resolution = (1280, 960), fontsize=18)

	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[2, 2] = GridLayout()
	gd = f[3, 1] = GridLayout()
	ge = f[3, 2] = GridLayout()

	ax = [
		Axis(ga[1,1], 
			xticks=(1:length(S), S_labels), 
			ylabel="CBI"),
		Axis(gb[1,1],
			yticks=0.0:25.0:100.0, 
			xlabel="Cue", 
			xticks=(1:length(C), C_labels), 
			ylabel="HA responses [%]"),
		Axis(gc[1,1], 
			xlabel="Cue", 
			xticks=(1:length(C), C_labels), 
			ylabel="Response time [sec]"),
		Axis(gd[1,1], 
			xlabel="Cue", 
			xticks=(1:length(C), C_labels), 
			ylabel="Omissions [%]"),
		Axis(ge[1,1], 
			xticks=(1:length(S), S_labels), 
			ylabel="Premature responses [%]")
		]

	plot_CBI!(ax[1], df, S)
	plot_responses!(ax[2], df, 2, C, S, S_labels)
	plot_RT!(ax[3], df, C, S, S_labels)
	plot_omissions!(ax[4], df, C, S, S_labels)
	plot_prematures!(ax[5], df, S)

	xlims!(ax[1], 0.5, length(S)+0.5)
	xlims!(ax[5], 0.5, length(S)+0.5)
	ylims!(ax[1], -1.4, 1.4)

	f[1, 2] = Legend(f, ax[2], framevisible=false, unique=true, tellwidth=false, tellheight=false)

	for (label, layout) in zip(["A", "B", "C", "D", "E"], [ga, gb, gc, gd, ge])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 26,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	save(string("./figures/", batch_ID, "_", study, ".eps"),f)
end

figure_probe("HO2", 
			"1vs1", 
			["1vs1_2"],
			["Second session"], 
			["2 kHz", "5 kHz", "8 kHz"]; 
			ID_excluded=["HO2_4", "HO2_11"], 
			S_excluded=["1vs1_1"]
			)

figure_probe("HO2", 
			"4vs1", 
			["4vs1_1", "4vs1_2"], 
			["First session", "Second session"], 
			["HC", "AC", "LC"]
			)

figure_drug("HO2", 
			"24hket",
			["Veh", "Ket_24h"],
			["Vehicle", "Ketamine 24h"], 
			["HC", "AC", "LC"]
			)

figure_drug("HO2", 
			"amph", 
			["Veh", "Amph"], 
			["Vehicle", "Amphetamine"], 
			["HC", "AC", "LC"]
			)

figure_drug("HO2", 
			"ket", 
			["Veh", "Ket"], 
			["Vehicle", "Ketamine"], 
			["HC", "AC", "LC"]
			)

figure_drug("HO2",
			"ket_ds", 
			["Veh", "Ket_03", "Ket_1", "Ket_3"], 
			["Vehicle", "Ketamine \n 0.3 mg/kg", "Ketamine \n 1 mg/kg", "Ketamine \n 3 mg/kg"], 
			["HC", "AC", "LC"]
			)

figure_probe("HO3", 
			"1vs1", 
			["1vs1_1", "1vs1_2"],
			["First session", "Second session"], 
			["2 kHz \n + lever lights", "lever lights", "8 kHz \n + lever lights"]
			)

figure_drug("HO3", 
			"ket_1vs1", 
			["Veh", "Ket"], 
			["Vehicle", "Ketamine"], 
			["2 kHz \n + lever lights", "lever lights", "8 kHz \n + lever lights"]
			)

figure_probe("HO3", 
			"1vs1_reversed", 
			["1vs1_reversed_1", "1vs1_reversed_2"],
			["First session", "Second session"], 
			["2 kHz \n + lever lights", "lever lights", "8 kHz \n lever lights"];
			reversed=true
			)