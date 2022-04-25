### A Pluto.jl notebook ###
# v0.19.0

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

Given how densely interconnected assemblies are, another expected property is that **activating a subset will result in activating the entire pattern**.
We will activate a random, small subset of the assembly and note how many other neurons of the assembly become active on the next timestep due to the recurrent connections between them (no learning).

First, we create a region and generate an assembly `m`.
"""

# ╔═╡ 6af41497-63d1-40b5-a8ad-f5c86b505789
begin
	A= Region(lib.params_input.sp, lib.params_input.tm);
	a= @chain bitrand(lib.Nin) A(_).active;
end

# ╔═╡ c068cf25-049a-4214-abfb-02268bc8000c
begin
	M= Region(lib.params_M.sp, lib.params_M.tm);
	assembly= lib.project!(M,a)
end

# ╔═╡ dfff8bd5-1399-43b8-9741-bd72f8042090
md"To begin with, we define a way to produce the partial stimulation:"

# ╔═╡ b1d3f0a8-8f9a-42b3-8019-5a071f5d682a
"`subset(x,p)` activates a random subset of neurons with ratio `p` from activation `x`. Note: `0≤p≤1`"
subset(x,p)= @chain [i for i= HierarchicalTemporalMemory.Truesof(x)] begin
	shuffle(_)[1 : round(Int,p*count(x))]
	HierarchicalTemporalMemory.bitarray(_,length(x))
end

# ╔═╡ 5cd6ac00-0de9-4266-ac87-4670e7d180a9
md"Note that, for every ratio:"

# ╔═╡ 6b065f0d-cacb-4590-be66-1ec7a6e6fc8c
count(subset(a,0.3)) == 0.3 * count(a)

# ╔═╡ 87caf94b-8b2e-4367-9795-2475163596d6
md"""
### Minicolumns

The `assembly` relies on the activation of local circuits by `a` through the strong feedforward connections.
These local circuits are called "minicolumns".
Deviating from the original assembly calculus, the assembly is strictly constrained to be a subpopulation of the active minicolumns, in general 1 neuron per minicolumn that wins a local competition.
The "votes" for this competition come from the region's recurrent connections.

To examine pattern completion we need the assembly's active minicolumns and also a subset of them:
"""

# ╔═╡ c17eb5f6-ae12-4907-a00c-ef01cfe854fb
begin
	c= M.sp(a)			# just the active minicolumns
	subsampling= 0.3
	ĉ= subset(c,subsampling)
end

# ╔═╡ 1e43ead6-f7fc-42cd-9911-c6d621ecd8f4
md"""
Note: the full activation of the region is
``$ M(a) = M.tm(M.sp(a)) = M.tm(c) $``
where `sp` means "spatial pooler" and `tm` means "temporal memory", the basic algorithms of HTM.
"""

# ╔═╡ 52b5365c-2c28-46fb-b1e7-0d5401e99ea0
md"""
### Partial stimulation

Every stimulation of a region will activate about the same number of neurons in total,
so what we mean by "partial stimulation" is that it should share a percentage of the full stimulation's neurons.

Since the full stimulation has cardinality (number of active neurons) $(count(M(a).active)) and, based on the previous rule, the subset should share $(subsampling * count(M(a).active)) neurons with it, we expect the overlap of the region's activation for each input to be close to the same value:
"""

# ╔═╡ 03f440f0-463e-4913-8c94-81285ce42d6f
lib.overlap( M.tm(c).active, M.tm(ĉ).active )

# ╔═╡ 466ebe2e-5ea4-4078-9264-db644ec18363
md"""
When the region is stimulated with the partial stimulation, we expect the recurrent connections within the layer to stimulate other neurons of the same assembly due to the interconnection property that was shown in [01-projection.jl](./open?path=notebooks/01-intro_projection.jl).

However, HTM dynamics won't allow those neurons to actually activate, because their minicolumns aren't active.
Instead, they will enter a "predictive" state.
In it, they are ready to win the minicolumn competition if the feedforward stimulation (proximal synapses) happens to activate their minicolumn.
We can compare the predictive neurons with the assembly and expect a high overlap.

## Pattern completion

To understand if and how well the $(subsampling) subset of the activation can recover the rest we'll rely on the overlap of the predicted neurons with the assembly.
To evaluate the overlap, we're going to compare:

- overlap of the assembly with neurons predicted after subset activation
- overlap of the assembly with neurons predicted after full activation

to answer the following questions:

- The coverage of the assembly. How much can be recovered?
- The percentage of correctly predicted neurons. Were there mispredictions or insufficient predictions?
"""

# ╔═╡ 0d6ca41f-a5a0-43ea-9ba3-1eda4c2f6ee4
subpredict_assembly= lib.overlap( assembly, M.tm(ĉ).predictive )

# ╔═╡ a4ed4028-9cfc-4a0b-92ea-b04bdc0ac058
fullpredict_assembly= lib.overlap( assembly, M.tm(c).predictive )

# ╔═╡ 54186c54-cca4-4033-b606-6fcec6de470e
coverage_full= fullpredict_assembly / count(assembly) *100

# ╔═╡ d0541d26-dc92-4987-9e62-66a05289acb2
coverage_subset= subpredict_assembly / count(assembly) *100

# ╔═╡ 7ae2283e-5441-44d6-9950-04d42e28e781
misprediction= 1 - subpredict_assembly / count(M.tm(ĉ).predictive)

# ╔═╡ 81243c19-dc2b-49af-8f7d-5969dd52291d
md"""
When the assembly is stimulated, the predictive neurons it generates predict $(round(coverage_full))% of the assembly.
This is a hallmark of a converged assembly: it can fully predict its own activation.

Compared to this, the $(subsampling) subset predicts $(round(coverage_subset))% of the assembly.
In other words, it retains $(round(coverage_subset/coverage_full*100))% of the full activation's predictive power, allowing a × $(round(coverage_subset/coverage_full / subsampling, sigdigits=2)) higher subsampling of the active neurons before prediction deteriorates.

### Assembly recall as a function of subsampling

"""

# ╔═╡ 32e81cbb-7da0-4767-a616-a99e631863c5
begin
	subset_fractions= 0:.05:1
	ĉ_sweep= map(p-> subset(c,p), subset_fractions)
	subpredict_coverage_sweep= map(c_sub-> lib.overlap( assembly, M.tm(c_sub).predictive )/fullpredict_assembly*100, ĉ_sweep)
end

# ╔═╡ 98487e5e-dec9-4419-90b8-f258a96f5729
plot(subset_fractions, subpredict_coverage_sweep,
	ylabel= "assembly recall %", xlabel= "assembly subsampling", label=:none,
	title= "Assembly recall as a function of subsampling")

# ╔═╡ 801c61f2-8981-4875-b714-c87ce84a287d
md"""
The graph shows that the coverage/recall is very sensitive to the assembly's subsampling.
Essentially, the assembly is robust up to about 30% of its neurons being activated; the rest are redundant, and are going to be activated by alternative connections on the next timestep.

*This function of assembly calculus is preserved by the HTM model.*

#### Robustness as a function of projection convergence time?

TODO:
- Manually converge the assembly and show this plot every 5 convergence steps.
"""

# ╔═╡ Cell order:
# ╠═351fd01c-8047-4d61-89fb-ab05e3bea911
# ╠═71fc74c6-317a-4418-8c0b-2ade268f7e6c
# ╟─0eb3cfdc-0f7b-4e0b-aace-c0a7a3669ead
# ╠═6af41497-63d1-40b5-a8ad-f5c86b505789
# ╠═c068cf25-049a-4214-abfb-02268bc8000c
# ╟─dfff8bd5-1399-43b8-9741-bd72f8042090
# ╠═b1d3f0a8-8f9a-42b3-8019-5a071f5d682a
# ╟─5cd6ac00-0de9-4266-ac87-4670e7d180a9
# ╠═6b065f0d-cacb-4590-be66-1ec7a6e6fc8c
# ╟─87caf94b-8b2e-4367-9795-2475163596d6
# ╠═c17eb5f6-ae12-4907-a00c-ef01cfe854fb
# ╟─1e43ead6-f7fc-42cd-9911-c6d621ecd8f4
# ╟─52b5365c-2c28-46fb-b1e7-0d5401e99ea0
# ╠═03f440f0-463e-4913-8c94-81285ce42d6f
# ╟─466ebe2e-5ea4-4078-9264-db644ec18363
# ╠═0d6ca41f-a5a0-43ea-9ba3-1eda4c2f6ee4
# ╠═a4ed4028-9cfc-4a0b-92ea-b04bdc0ac058
# ╠═54186c54-cca4-4033-b606-6fcec6de470e
# ╠═d0541d26-dc92-4987-9e62-66a05289acb2
# ╠═7ae2283e-5441-44d6-9950-04d42e28e781
# ╟─81243c19-dc2b-49af-8f7d-5969dd52291d
# ╠═32e81cbb-7da0-4767-a616-a99e631863c5
# ╠═98487e5e-dec9-4419-90b8-f258a96f5729
# ╟─801c61f2-8981-4875-b714-c87ce84a287d
