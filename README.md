# JBT

This repository contains code used for each part of the analysis pipeline for data coming from the Judgement Bias Task (JBT). The JBT is an operant, decision-making task where subjects are asked to interpret an ambiguous stimulus and take an action according to its interpreted meaning. One interpretation could lead to a large reward, whereas the alternative leads to a much smaller reward. The task hypothesis is that the subjects' prior expectations of reward would bias their interpretation of the ambiguous cue and consequently lead to biased actions.

In the /src folder, there exist a `JBT` Module, written in Julia (v1.7) that contains code to read raw data, coming from the KLimbic software that runs the operant chambers, convert the data to DataFrames, plot useful summaries (such as accuracy, response times, an index of the subjects' interpretation bias etc), model the data using a hierarchical GLM statistical model, which takes into account the effect of previous trials on the currect actions and finally visualise the model results, in the form of posterior distribution plots.

The analysis scripts `model_fit.jl`, `model_comparison.jl`, `figures_analysis.jl` and `figures_exp.jl`, are examples of using the `JBT` Module to perform model fitting, model comparison and plotting of experimental and model outcomes. These constituted the results of 2 chapters of my PhD thesis.

In order to use the `JBT` Module, please clone this repository, go into the downloaded JBT folder 
```shell
cd JBT
```
and then open a `julia` instance and run 
```julia
import Pkg
Pkg.activate(".")
```

This will activate the `JBT` environment, along with all its dependencies.