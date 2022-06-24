using JBT
using Serialization
using ArviZ

include("utils.jl")

df = get_batch_df(
    Dict(
        "./exp/probe/baseline/baseline_effect/"=>"baseline",
        "./exp/probe/ketamine/ketamine_effect/"=>"ketamine",
        "./exp/probe/ketamine/vehicle_effect/"=>"vehicle"
        )
    )

chain_RR = deserialize("./chains/chain_RR_effect_batch.jls")
chain_learning = deserialize("./chains/chain_learning_effect_batch.jls")

idt_learning = chain_to_idt(df, chain_learning, cue_reward_rate)
idt_RR = chain_to_idt(df, chain_RR, reward_rate)

d_model = Dict("learning"=>idt_learning, "RR"=>idt_RR)
compare(d_model, ic="loo", method="stacking", scale="log")
