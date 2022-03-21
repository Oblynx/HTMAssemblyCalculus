### A Pluto.jl notebook ###
# v0.18.2

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

# ╔═╡ ae5a3aa6-4071-4ab2-8130-700b583ec60d
md"""
# Association of 2 assemblies

Using the **projection** operation shown in [01-intro_projection.jl](./open?path=notebooks/01-intro_projection.jl), we are now going to create 2 assemblies in the same region and see what happens if they fire *at the same time*.

### Regions

We'll use 3 regions in this experiment:

- 2 "input" regions A,B. These identical regions have externally activated sets of neurons.
- 1 "merge" region M. This is where the 2 assemblies will be generated and interact.

## Generating the assemblies

We'll copy the Region configuration of [01-intro_projection.jl](./open?path=notebooks/01-intro_projection.jl):
"""

# ╔═╡ e3bda23f-6ec8-4977-9fad-e7cf99aa8bb3
md"The 3 regions:"

# ╔═╡ d5286246-3265-421e-b724-13b963f09887
A= Region(lib.params_input.sp, lib.params_input.tm);

# ╔═╡ 1420b89c-4bd5-4278-8a70-92304a9a5723
B= Region(lib.params_input.sp, lib.params_input.tm);

# ╔═╡ e729ea16-387f-439d-843a-5383ad93703e
M= Region(lib.params_M.sp, lib.params_M.tm);

# ╔═╡ 5e7753e3-794a-4aa1-ae4a-a9cbedc9fc52
md"On regions A,B we can statically generate the input neuron activations `a,b`. We'll use `a,b` to stimulate region M and generate assemblies there. Input regions don't have to adapt."

# ╔═╡ eb32b26e-7b2e-4d0e-a5fa-81484180f07e
approxSubsample!(x,p)= begin
	idx= randsubseq(findall(x), 1-p)
	x[idx].= false
	x
end

# ╔═╡ fef94d49-27e4-4809-9804-270e9270f958
a= @chain bitrand(lib.Nin) A(_).active approxSubsample!(_, 1/4);

# ╔═╡ 5fff0d6b-5098-4500-a521-c0553ff5a4ef
b= @chain bitrand(lib.Nin) B(_).active approxSubsample!(_, 1/4);

# ╔═╡ 52f798a4-a12c-462e-b07f-80a15cf7d159
md"We can now project them to M and generate the assemblies `â,b̂`:"

# ╔═╡ 0b396e36-3022-4c39-8e44-84ebd42d7b19
begin
	reset!(M)
	â₀= lib.project!(M,a)
	b̂₀= lib.project!(M,b)
end

# ╔═╡ 9f0f1131-52cc-4894-a1e0-7999fbd5234a
md"""
## Assembly overlap

The sparsity of neuronal activations ensures that only a very small proportion of a region's neurons will take part in each assembly.
This implies that 2 unrelated assemblies have very low probability to overlap even at a single neuron.
Conversely, a high overlap is a strong indication that 2 assemblies are not unrelated.
"""

# ╔═╡ 0e31f16f-f2c3-4497-ad0b-47d0fc393d4c
"`overlap(x,y)` of 2 assemblies in the same region is the number of neurons that are active in both. Symbol: `x&y`"
overlap(x,y)= count(x .& y)

# ╔═╡ b2e2c302-eefd-41ec-97ca-aa302d62f2f1
md"""
Right after creating the assemblies, their overlap is $(overlap(â₀,b̂₀)) as expected.
What will happen if the region receives both inputs `a,b` at the same time consistently?

Given the sparsity of neuron activations, the assemblies `â,b̂` will be competing for the same few active neurons.
Let the resulting activation be `âb̂`. We would expect it to be a "compromise" between `â,b̂`, overlapping by similar amounts with each of them.
As `âb̂` becomes an assembly after repeated activations, do the assemblies `â` and `b̂` remain unaffected?

To answer this question we will repeatedly stimulate `M` with `a|b` (`OR` of a,b), producing `âb̂ᵢ`.
After each stimulation, we will also stimulate (without learning) with `a` and `b` separately, in order to measure:

- the overlap of `âb̂ᵢ` with the initial assemblies `â₀` and `b̂₀`. This is our reference frame: whether the activation remains balanced.
- the overlap of `âᵢ` and `b̂ᵢ` with `â₀` and `b̂₀` respectively. This will show whether the original assemblies are changing.
- the overlap of `âᵢ` and `b̂ᵢ`

The latter is the study focus. If this overlap increases, we will say that `â` and `b̂` are **associating**.
"""

# ╔═╡ 2cc30337-2c48-476f-b45f-c2bfb788a3db
measure(M, âb̂)= begin
	â= lib.probe(M,a).active
	b̂= lib.probe(M,b).active
	(
		ref= (overlap(âb̂,â₀), overlap(âb̂,b̂₀)),
		Δa= (overlap(â,â₀), overlap(b̂,b̂₀)),
		assoc= overlap(â,b̂)
	)
end

# ╔═╡ 8330fc10-22c9-419b-887b-7d91ee580dc3
stimulateMeasure!(M)= begin
	# reset predictive state to activation from a.|b
	# TODO: step!(R,x, learn= false)
	#step!(M, a.|b, falses(0), false)
	âb̂= step!(M, a.|b).active
	measure(M, âb̂)
end

# ╔═╡ 76758fe0-0136-44dc-9144-193aea2130f5
md"With the measurements defined, we can now run the experiment."

# ╔═╡ 39437dfc-38cb-412f-8fc6-02fa3f646473
md"""
## Results

The figure below shows how the assemblies are diverging from their original set until $t=6$, while at the same they are approaching closer together:
"""

# ╔═╡ 47c839f0-d519-4267-ac2b-535338fae0c4
T= 50   # time to converge

# ╔═╡ 13e856d8-2eaf-4326-8026-b4e7edebf4c0
begin
	M_assoc= deepcopy(M)
	Random.seed!(0)
	measurements= [
		measure(M_assoc, M_assoc(a.|b).active);    # starting point
		[stimulateMeasure!(M_assoc) for t= 1:T]
	]
end

# ╔═╡ 3660c5f0-cc2a-4b95-86b5-2392f5bc16ec
â= lib.probe(M_assoc,a).active;

# ╔═╡ a9d8b7a7-352e-4b65-931a-006958e9d465
b̂= lib.probe(M_assoc,b).active;

# ╔═╡ d109633b-1740-4da4-ad10-eb2f2f6527de
begin
	p= get_color_palette(:auto, plot_color(:white))
	p2= vec([0.8,1]*p[1:3]')
	plot(map(f-> f.(measurements), [x-> x.Δa[1], x-> x.Δa[2], x-> x.ref[1], x-> x.ref[2], x-> x.assoc]),
		linecolor= p2[1:5]', linewidth= 1.6,
		title="Association:\n changing co-occurrent assemblies",
		xlabel="t", ylabel="ovelap", legend=:right,
		labels= hcat("âᵢ & â₀","b̂ᵢ & b̂₀",
			"âb̂ᵢ & â₀", "âb̂ᵢ & b̂₀",
			"âᵢ & b̂ᵢ"),
		minorgrid= true,
	)
end

# ╔═╡ cd469ae2-74d4-439e-ad27-9f2ea999a13c
md"""
The overlap $\texttt{âᵢ} \& \texttt{b̂ᵢ}$ is at least **$(
round(minimum(map(x->x.assoc, measurements[end-10:end])) / mean([count(â),count(b̂)]) * 100)|> Int
)%** of the number of active neurons.
This is comparable to the figure reported in the paper (8-10%).

The blue, cyan lines show that the original assemblies are changing due to the co-occurrence.

Furthermore, despite its evolution, `âb̂` remains balanced between the 2 "source" assemblies, as shown by the similarity of the red and orange lines.

#### Co-occurrence → co-adaptation

The change we observed in co-occurrent assemblies is called _association_.

The associated assemblies in M are `â`, `b̂`.
"""

# ╔═╡ cf9cc6f2-e74f-4bf3-b936-2aa436c32f33
md"Their overlap is $(round(overlap(â,b̂)/mean([count(â),count(b̂)]) * 100)|> Int
)%."

# ╔═╡ e3583c97-962b-4517-b5ca-b8f469fab5ec
md"""
## Is association conserved through projection?

Let's project `a,b` to a new region and see if the resulting assemblies have the same association, without any extra training.
In such a case we would say that the association between the assemblies is conserved when they are viewed in a different brain area.

The new region:
"""

# ╔═╡ 07c64927-66bd-43d9-b1c9-d05414ceb8df
P= Region(lib.params_M.sp, lib.params_M.tm);

# ╔═╡ aacd1964-b120-4bb4-8be5-a4b1f2908f82
begin
	reset!(P)
	a_p= lib.project!(P,â)
	b_p= lib.project!(P,b̂)
end

# ╔═╡ 54753e1d-fa82-4bab-9ce5-f1c81cc3e133
md"""
The overlap of the resulting assemblies is **$(overlap(a_p,b_p))**.

This property of assembly calculus is **not satisfied** by this HTM model.
"""

# ╔═╡ 8272eb47-1d34-4ab0-81b9-9e558555c5ca
md"""
### Why is the association not conserved?

Since `â`, `b̂` have significant overlap, it is likely that their projections have some minicolumns in common.
With the following functions we can calculate how many minicolumns are activated in common.
"""

# ╔═╡ 31e6d515-bbe7-45e3-ac29-14bec5667a73
"`minicolumnOverlap(a,b,R)` is the number of minicolumns that are active in both a,b (of region R)"
minicolumnOverlap(a,b,R)= overlap(
	any(reshape(a,R.tm.params.k,:),dims=1),
	any(reshape(b,R.tm.params.k,:),dims=1)
)

# ╔═╡ a735fd78-186a-442b-bbc1-784b296c05b8
"`minicolumnOverlapFraction(a,b,R)` normalizes the number of overlapping minicolumns by the number of minicolumns expected to fire in each activation"
minicolumnOverlapFraction(a,b,R)= minicolumnOverlap(a,b,R) / (R.tm.params.Nc*R.sp.params.s)

# ╔═╡ b31a588e-9b36-4826-9f8f-0277af0159ee
md"If we activate the region `P` with `â` and `b̂` (before their projection), we get only a very small fraction of overlapping minicolumns."

# ╔═╡ b207362f-575d-4d0f-83f1-1b9100739a21
begin
	reset!(P)
	a_project_p= lib.probe(P,â).active
	b_project_p= lib.probe(P,b̂).active
end

# ╔═╡ b3220822-1b9c-4b81-b2a0-2673f0c4fc8f
md"Fraction of overlapping minicolumns: **$(round(minicolumnOverlapFraction(a_project_p, b_project_p, P)*100, digits=1))%**."

# ╔═╡ 90229641-1812-458f-9517-d8898b0eaf31
md"If we investigate the assemblies, we see that the assemblies `â, b̂` are represented by very different counts of neurons, meaning that `b` was never learned as strongly. This imbalance might be the cause of the association not being conserved."

# ╔═╡ 20a84dba-21ff-4222-a832-eb4c54e38638
count(a), count(b)

# ╔═╡ 8e8add95-758b-4984-ad51-e5516bb36c66
count(â), count(b̂)

# ╔═╡ 785a5a7c-78f2-4700-ae5c-60413f129219
count(a_p), count(b_p)

# ╔═╡ Cell order:
# ╠═351fd01c-8047-4d61-89fb-ab05e3bea911
# ╠═71fc74c6-317a-4418-8c0b-2ade268f7e6c
# ╟─ae5a3aa6-4071-4ab2-8130-700b583ec60d
# ╟─e3bda23f-6ec8-4977-9fad-e7cf99aa8bb3
# ╠═d5286246-3265-421e-b724-13b963f09887
# ╠═1420b89c-4bd5-4278-8a70-92304a9a5723
# ╠═e729ea16-387f-439d-843a-5383ad93703e
# ╟─5e7753e3-794a-4aa1-ae4a-a9cbedc9fc52
# ╟─eb32b26e-7b2e-4d0e-a5fa-81484180f07e
# ╠═fef94d49-27e4-4809-9804-270e9270f958
# ╠═5fff0d6b-5098-4500-a521-c0553ff5a4ef
# ╟─52f798a4-a12c-462e-b07f-80a15cf7d159
# ╠═0b396e36-3022-4c39-8e44-84ebd42d7b19
# ╟─9f0f1131-52cc-4894-a1e0-7999fbd5234a
# ╟─0e31f16f-f2c3-4497-ad0b-47d0fc393d4c
# ╟─b2e2c302-eefd-41ec-97ca-aa302d62f2f1
# ╠═2cc30337-2c48-476f-b45f-c2bfb788a3db
# ╠═8330fc10-22c9-419b-887b-7d91ee580dc3
# ╟─76758fe0-0136-44dc-9144-193aea2130f5
# ╟─39437dfc-38cb-412f-8fc6-02fa3f646473
# ╠═47c839f0-d519-4267-ac2b-535338fae0c4
# ╠═13e856d8-2eaf-4326-8026-b4e7edebf4c0
# ╠═3660c5f0-cc2a-4b95-86b5-2392f5bc16ec
# ╠═a9d8b7a7-352e-4b65-931a-006958e9d465
# ╟─d109633b-1740-4da4-ad10-eb2f2f6527de
# ╟─cd469ae2-74d4-439e-ad27-9f2ea999a13c
# ╠═cf9cc6f2-e74f-4bf3-b936-2aa436c32f33
# ╟─e3583c97-962b-4517-b5ca-b8f469fab5ec
# ╠═07c64927-66bd-43d9-b1c9-d05414ceb8df
# ╠═aacd1964-b120-4bb4-8be5-a4b1f2908f82
# ╟─54753e1d-fa82-4bab-9ce5-f1c81cc3e133
# ╟─8272eb47-1d34-4ab0-81b9-9e558555c5ca
# ╟─31e6d515-bbe7-45e3-ac29-14bec5667a73
# ╟─a735fd78-186a-442b-bbc1-784b296c05b8
# ╠═b31a588e-9b36-4826-9f8f-0277af0159ee
# ╠═b207362f-575d-4d0f-83f1-1b9100739a21
# ╟─b3220822-1b9c-4b81-b2a0-2673f0c4fc8f
# ╟─90229641-1812-458f-9517-d8898b0eaf31
# ╠═20a84dba-21ff-4222-a832-eb4c54e38638
# ╠═8e8add95-758b-4984-ad51-e5516bb36c66
# ╠═785a5a7c-78f2-4700-ae5c-60413f129219
