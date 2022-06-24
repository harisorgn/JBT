
@model function logistic_past(choices, data::JBTExperiment, ::Type{T} = Float64) where {T}

	if choices === missing
		choices = Vector{T}(undef, length(data.ID))
	end

	μ_b ~ filldist(Normal(0, 1), data.n_conditions)
	σ_b ~ filldist(Exponential(1), data.n_conditions)
	b_norm ~ filldist(Normal(0, 1), data.n_conditions, data.n_subjects)
	b = @. μ_b + b_norm * σ_b

	μ_r ~ filldist(Normal(0, 1), data.n_conditions)
	σ_r ~ filldist(Exponential(1), data.n_conditions)
	w_r_norm ~ filldist(Normal(0, 1), data.n_conditions, data.n_subjects)
	w_r = @. μ_r + w_r_norm * σ_r

	μ_p ~ filldist(Normal(0, 1), data.n_conditions)
	σ_p ~ filldist(Exponential(1), data.n_conditions)
	w_p_norm ~ filldist(Normal(0, 1), data.n_conditions, data.n_subjects)
	w_p = @. μ_p + w_p_norm * σ_p

	μ_ε ~ Normal(-1.5, 0.2)
	σ_ε ~ Exponential(0.5)
	ε_norm ~ filldist(Normal(0,1), data.n_subjects)
	ε = @. cdf(Normal(0,1), μ_ε + ε_norm * σ_ε)

	A = map((ID, condition, I_m1, R) ->  b[condition,ID] + w_p[condition,ID]*I_m1 + w_r[condition,ID]*R,
			data.ID, data.conditions, data.I_m1, data.R)

	P = @. (1.0-ε[data.ID])*(1.0/(1.0 + exp(-A))) + ε[data.ID]/2.0

	for i in eachindex(choices)
		choices[i] ~ Binomial(1, P[i])
	end

	return choices
end
