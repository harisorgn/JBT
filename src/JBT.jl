module JBT

export read_JBT, df_to_JBT, JBTExperiment,
        plot_CBI!, plot_responses!, plot_RT!, plot_omissions!, plot_prematures!,
        plot_diff_CBI!, condition_dict, logistic_past, chain_to_idt

using DataFrames
using CSV
using CairoMakie
using ColorSchemes
using LaTeXStrings
using Distributions: Uniform
using HypothesisTests
using Serialization
using ArviZ: hdi
using KernelDensity: kde
using Statistics: mean, std
using LinearAlgebra
using Turing

include("read.jl")
include("JBTExperiment.jl")
include("metrics.jl")
include("plot.jl")
include("model.jl")
include("model_compare.jl")

end
