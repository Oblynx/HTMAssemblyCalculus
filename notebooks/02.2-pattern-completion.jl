### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 351fd01c-8047-4d61-89fb-ab05e3bea911
begin
	import Pkg; Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using HierarchicalTemporalMemory
	using Random, Chain, Setfield, Statistics, Plots, PlutoUI
end

# ╔═╡ 71fc74c6-317a-4418-8c0b-2ade268f7e6c
"Include assembly operations library."
module lib
	include(joinpath(dirname(Base.active_project()), "src/assembly_operations.jl"))
end

# ╔═╡ 0eb3cfdc-0f7b-4e0b-aace-c0a7a3669ead
md"""
# Pattern completion

Given how densely interconnected assemblies are, another expected property is that they can complete their own activation.
We will activate a random, small subset of the assembly and note how many other neurons of the assembly become active on the next timestep due to the recurrent connections between them (no learning).

To begin with, we define a way to produce the partial stimulation:
"""

# ╔═╡ b1d3f0a8-8f9a-42b3-8019-5a071f5d682a
"`subset(x,p)` activates a random subset of neurons with ratio `p` from activation `x`. Note: `0≤p≤1`"
subset(x,p)= @chain [i for i= HierarchicalTemporalMemory.Truesof(x)] begin
	shuffle(_)[1 : Int(p*count(x))]
	HierarchicalTemporalMemory.bitarray(_,length(a))
end

# ╔═╡ 5cd6ac00-0de9-4266-ac87-4670e7d180a9
md"Note that:"

# ╔═╡ 6b065f0d-cacb-4590-be66-1ec7a6e6fc8c
count(subset(a,0.2)) == 0.2 * count(a)

# ╔═╡ 87caf94b-8b2e-4367-9795-2475163596d6
md"""
From the input `a` we can now produce the full stimulation of the region `M(a)` and also the partial stimulation `M(subset(a,0.2))`.
Every stimulation of a region will activate about the same number of neurons in total,
so what we mean by "partial stimulation" is that it should share a percentage of the full stimulation's neurons.
Since the full stimulation has cardinality (number of active neurons) $(count(M(a).active)) and, based on the previous rule, the 0.2 subset should share $(0.2count(M(a).active)) neurons with it, we expect the following overlap to be close to the same value:
"""

# ╔═╡ 03f440f0-463e-4913-8c94-81285ce42d6f
overlap( M(a).active, M(subset(a,0.2)).active )

# ╔═╡ 2253579f-69d9-47ea-bffb-f06a63785222
md"""

"""

# ╔═╡ 37b7c578-83ac-4d74-bda4-e6cafd0ad000
M.tm(M.sp(a), M(a).predictive)

# ╔═╡ 8ce48b65-09a0-4375-942d-32459ae79cd9
begin
    _M2= deepcopy(M)
	x= subset(a,0.2)
	α= [step!(_M2, x).active for t=1:30]
	plot( overlap.( Ref(_M2(a).active), α) )
end

# ╔═╡ Cell order:
# ╠═351fd01c-8047-4d61-89fb-ab05e3bea911
# ╟─71fc74c6-317a-4418-8c0b-2ade268f7e6c
# ╠═0eb3cfdc-0f7b-4e0b-aace-c0a7a3669ead
# ╟─b1d3f0a8-8f9a-42b3-8019-5a071f5d682a
# ╟─5cd6ac00-0de9-4266-ac87-4670e7d180a9
# ╠═6b065f0d-cacb-4590-be66-1ec7a6e6fc8c
# ╟─87caf94b-8b2e-4367-9795-2475163596d6
# ╠═03f440f0-463e-4913-8c94-81285ce42d6f
# ╠═2253579f-69d9-47ea-bffb-f06a63785222
# ╠═37b7c578-83ac-4d74-bda4-e6cafd0ad000
# ╠═8ce48b65-09a0-4375-942d-32459ae79cd9
