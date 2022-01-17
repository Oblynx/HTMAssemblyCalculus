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
# ╟─ae5a3aa6-4071-4ab2-8130-700b583ec60d
# ╠═77db44b5-01bd-46c3-afa5-60fca3743d92
# ╟─e3bda23f-6ec8-4977-9fad-e7cf99aa8bb3
# ╠═d5286246-3265-421e-b724-13b963f09887
# ╠═1420b89c-4bd5-4278-8a70-92304a9a5723
# ╠═e729ea16-387f-439d-843a-5383ad93703e
# ╟─5e7753e3-794a-4aa1-ae4a-a9cbedc9fc52
# ╠═fef94d49-27e4-4809-9804-270e9270f958
# ╠═5fff0d6b-5098-4500-a521-c0553ff5a4ef
# ╟─52f798a4-a12c-462e-b07f-80a15cf7d159
# ╠═0b396e36-3022-4c39-8e44-84ebd42d7b19
# ╟─9f0f1131-52cc-4894-a1e0-7999fbd5234a
# ╟─0e31f16f-f2c3-4497-ad0b-47d0fc393d4c
# ╟─b2e2c302-eefd-41ec-97ca-aa302d62f2f1
# ╠═2cc30337-2c48-476f-b45f-c2bfb788a3db
# ╟─8330fc10-22c9-419b-887b-7d91ee580dc3
# ╟─76758fe0-0136-44dc-9144-193aea2130f5
# ╟─39437dfc-38cb-412f-8fc6-02fa3f646473
# ╠═13e856d8-2eaf-4326-8026-b4e7edebf4c0
# ╟─d109633b-1740-4da4-ad10-eb2f2f6527de
# ╟─cd469ae2-74d4-439e-ad27-9f2ea999a13c
# ╠═0eb3cfdc-0f7b-4e0b-aace-c0a7a3669ead
# ╟─b1d3f0a8-8f9a-42b3-8019-5a071f5d682a
# ╟─5cd6ac00-0de9-4266-ac87-4670e7d180a9
# ╠═6b065f0d-cacb-4590-be66-1ec7a6e6fc8c
# ╟─87caf94b-8b2e-4367-9795-2475163596d6
# ╠═03f440f0-463e-4913-8c94-81285ce42d6f
# ╠═2253579f-69d9-47ea-bffb-f06a63785222
# ╠═37b7c578-83ac-4d74-bda4-e6cafd0ad000
# ╠═8ce48b65-09a0-4375-942d-32459ae79cd9
