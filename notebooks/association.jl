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

# ╔═╡ ae5a3aa6-4071-4ab2-8130-700b583ec60d
md"""
# Association of 2 assemblies

Using the **projection** operation shown in [projection.jl](./open?path=notebooks/projection.jl), we are now going to create 2 assemblies in the same region and see what happens if they fire *at the same time*.

### Regions

We'll use 3 regions in this experiment:

- 2 "input" regions A,B. These identical regions have externally activated sets of neurons.
- 1 "merge" region M. This is where the 2 assemblies will be generated and interact.

## Generating the assemblies

We'll copy the Region configuration of [projection.jl](./open?path=notebooks/projection.jl):
"""

# ╔═╡ 77db44b5-01bd-46c3-afa5-60fca3743d92
begin
	Nin= 1e3|> Int         # input size
	Nn= 4e5                # number of neurons in each area
	k= 10                  # neurons per minicolumn
	T= 30                  # steps to convergence
	_Nc()= floor(Int,Nn/k) # number of minicolumns
	sparsity()= 1/sqrt(_Nc());
	params_input= (
		sp= SPParams(szᵢₙ=Nin, szₛₚ=_Nc(), s= sparsity(), prob_synapse=1e-3, enable_local_inhibit=false),
		tm= TMParams(Nc=_Nc(), k=k, p⁺_01= .11, θ_stimulus_learn=15, θ_stimulus_activate= 22)
	)
	params_M= @set params_input.sp.szᵢₙ= params_input.tm.Nₙ
end

# ╔═╡ e3bda23f-6ec8-4977-9fad-e7cf99aa8bb3
md"The 3 regions:"

# ╔═╡ d5286246-3265-421e-b724-13b963f09887
A= Region(params_input.sp, params_input.tm);

# ╔═╡ 1420b89c-4bd5-4278-8a70-92304a9a5723
B= Region(params_input.sp, params_input.tm);

# ╔═╡ e729ea16-387f-439d-843a-5383ad93703e
M= Region(params_M.sp, params_M.tm);

# ╔═╡ 5e7753e3-794a-4aa1-ae4a-a9cbedc9fc52
md"On regions A,B we can statically generate the input neuron activations `a,b`. We'll use `a,b` to stimulate region M and generate assemblies there. Input regions don't have to adapt."

# ╔═╡ fef94d49-27e4-4809-9804-270e9270f958
a= @chain bitrand(Nin) A(_).active;

# ╔═╡ 5fff0d6b-5098-4500-a521-c0553ff5a4ef
b= @chain bitrand(Nin) B(_).active;

# ╔═╡ 52f798a4-a12c-462e-b07f-80a15cf7d159
md"We can now project them to M and generate the assemblies `am,bm`:"

# ╔═╡ 0b396e36-3022-4c39-8e44-84ebd42d7b19
begin
	reset!(M)
    am₀= lib.project!(M,a)
	bm₀= lib.project!(M,b)
end

# ╔═╡ 9f0f1131-52cc-4894-a1e0-7999fbd5234a
md"""
## Assembly overlap

The sparsity of neuronal activations ensures that only a very small proportion of a region's neurons will take part in each assembly.
This implies that 2 unrelated assemblies have very low probability to overlap even at a single neuron.
Conversely, a high overlap is a strong indication that 2 assemblies are not unrelated.
"""

# ╔═╡ 0e31f16f-f2c3-4497-ad0b-47d0fc393d4c
"`overlap(x,y)` of 2 assemblies in the same region is the number of neurons that are active in both"
overlap(x,y)= count(x .& y)

# ╔═╡ b2e2c302-eefd-41ec-97ca-aa302d62f2f1
md"""
Right after creating the assemblies, their overlap is $(overlap(am₀,bm₀)) as expected.
What will happen if the region receives both inputs `a,b` at the same time consistently?

Given the sparsity of neuron activations, the assemblies `am,bm` will be competing for the same few active neurons.
Let the resulting activation be `abm`. We would expect it to be a "compromise" between `am,bm`, overlapping by similar amounts with each of them.
As `abm` becomes an assembly after repeated activations, do the assemblies `am` and `bm` remain unaffected?

To answer this question we will repeatedly stimulate `M` with `a|b` (`OR` of a,b), producing `abmᵢ`.
After each stimulation, we will also stimulate (without learning) with `a` and `b`, in order to measure:

- the overlap of `abmᵢ` with the initial assemblies `am₀` and `bm₀`. This is our reference frame: whether the activation remains balanced.
- the overlap of `amᵢ` and `bmᵢ` with `am₀` and `bm₀` respectively. This will show whether the original assemblies are changing.
- the overlap of `amᵢ` and `bmᵢ`

The latter is the study focus. If this overlap increases, we will say that `am` and `bm` are **associating**.
"""

# ╔═╡ 2cc30337-2c48-476f-b45f-c2bfb788a3db
measure(M, abm)= begin
	am= M(a).active
	bm= M(b).active
	(
		ref= (overlap(abm,am₀), overlap(abm,bm₀)),
		change= (overlap(am,am₀), overlap(bm,bm₀)),
		assoc= overlap(am,bm)
	)
end

# ╔═╡ 8330fc10-22c9-419b-887b-7d91ee580dc3
stimulateMeasure!(M)= begin
	abm= step!(M, a.|b).active
	measure(M, abm)
end

# ╔═╡ 13e856d8-2eaf-4326-8026-b4e7edebf4c0
begin
	_M= deepcopy(M)
	measurements= [stimulateMeasure!(_M) for t= 1:T]
	md"With the measurements defined, we can now run the experiment."
end

# ╔═╡ d109633b-1740-4da4-ad10-eb2f2f6527de
plot(map(f-> f.(measurements), [x-> x.change[1], x-> x.change[2], x-> x.assoc]))

# ╔═╡ 3097e22b-d4a2-42b5-90de-2b295d95cfd0


# ╔═╡ cd469ae2-74d4-439e-ad27-9f2ea999a13c
md"The overlap rises to about 10% as reported in the paper. This is the association between the 2 assemblies."

# ╔═╡ Cell order:
# ╠═351fd01c-8047-4d61-89fb-ab05e3bea911
# ╠═71fc74c6-317a-4418-8c0b-2ade268f7e6c
# ╟─ae5a3aa6-4071-4ab2-8130-700b583ec60d
# ╠═77db44b5-01bd-46c3-afa5-60fca3743d92
# ╟─e3bda23f-6ec8-4977-9fad-e7cf99aa8bb3
# ╠═d5286246-3265-421e-b724-13b963f09887
# ╠═1420b89c-4bd5-4278-8a70-92304a9a5723
# ╠═e729ea16-387f-439d-843a-5383ad93703e
# ╟─5e7753e3-794a-4aa1-ae4a-a9cbedc9fc52
# ╠═fef94d49-27e4-4809-9804-270e9270f958
# ╠═5fff0d6b-5098-4500-a521-c0553ff5a4ef
# ╠═52f798a4-a12c-462e-b07f-80a15cf7d159
# ╠═0b396e36-3022-4c39-8e44-84ebd42d7b19
# ╟─9f0f1131-52cc-4894-a1e0-7999fbd5234a
# ╠═0e31f16f-f2c3-4497-ad0b-47d0fc393d4c
# ╟─b2e2c302-eefd-41ec-97ca-aa302d62f2f1
# ╠═2cc30337-2c48-476f-b45f-c2bfb788a3db
# ╠═8330fc10-22c9-419b-887b-7d91ee580dc3
# ╠═13e856d8-2eaf-4326-8026-b4e7edebf4c0
# ╠═d109633b-1740-4da4-ad10-eb2f2f6527de
# ╠═3097e22b-d4a2-42b5-90de-2b295d95cfd0
# ╠═cd469ae2-74d4-439e-ad27-9f2ea999a13c
