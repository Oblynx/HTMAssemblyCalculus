### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 42619509-b547-4931-b0b3-82d02d78ad94
begin
	import Pkg; Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using HierarchicalTemporalMemory
	using Random, Chain, Setfield, Statistics, Plots, PlutoUI
	#using PlotlyJS;	plotlyjs()
end

# ╔═╡ 9c8aafe6-e978-47d3-893e-b3729b2b9439
md"""
# Neuron assemblies

> Assemblies are large populations of neurons believed to imprint memories, concepts, words, and other cognitive information.

This concept was introduced by \\cite

[1]: "Brain computation by assemblies of neurons" https://www.pnas.org/content/117/25/14464
"""

# ╔═╡ Cell order:
# ╠═42619509-b547-4931-b0b3-82d02d78ad94
# ╠═9c8aafe6-e978-47d3-893e-b3729b2b9439
