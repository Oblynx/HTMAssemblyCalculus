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
	#using PlotlyJS;	plotlyjs()
end

# ╔═╡ ae5a3aa6-4071-4ab2-8130-700b583ec60d
md"""
# Association of 2 assemblies

Using the **projection** operation shown in [projection.jl](./open?path=notebooks/projection.jl), we are now going to create 2 assemblies in the same region and see what happens if they fire *at the same time*.

### Regions

We'll use 3 regions in this experiment:

- 2 "input" regions A,B. These identical regions have externally activated sets of neurons.
- 1 "merge" region M. This is where the 2 assemblies will interact.

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
md"On regions A,B we can statically generate the input neuron activations `a,b`. We'll use `a,b` to stimulate region M and generate assemblies there."

# ╔═╡ fef94d49-27e4-4809-9804-270e9270f958
a= @chain bitrand(Nin) A(_).active;

# ╔═╡ 5fff0d6b-5098-4500-a521-c0553ff5a4ef
b= @chain bitrand(Nin) B(_).active;

# ╔═╡ f0bd22ad-8e33-457e-b9b5-82805a3cc547
project(x,r)= Iterators.drop((step!(r,x).active for t= 1:T), T-1)|> first

# ╔═╡ 0b396e36-3022-4c39-8e44-84ebd42d7b19
begin
    am= project(a,M)
	bm= project(b,M)
end

# ╔═╡ 9f0f1131-52cc-4894-a1e0-7999fbd5234a
md"## Overlap"

# ╔═╡ b2e2c302-eefd-41ec-97ca-aa302d62f2f1
md"Right after creating the assemblies, their overlap is $(count(am.&bm)). We'll cause both input activations to fire simultaneously and see the result on the overlap."

# ╔═╡ 2cc30337-2c48-476f-b45f-c2bfb788a3db
abm= project(a.|b, M)

# ╔═╡ e7bb6538-e885-4ff9-8246-d888ced1da0f
am2, bm2= M(a).active, M(b).active

# ╔═╡ 8330fc10-22c9-419b-887b-7d91ee580dc3
count(am2 .& bm2)

# ╔═╡ cd469ae2-74d4-439e-ad27-9f2ea999a13c
md"The overlap rises to about 10% as reported in the paper. This is the association between the 2 assemblies."

# ╔═╡ Cell order:
# ╠═351fd01c-8047-4d61-89fb-ab05e3bea911
# ╟─ae5a3aa6-4071-4ab2-8130-700b583ec60d
# ╠═77db44b5-01bd-46c3-afa5-60fca3743d92
# ╟─e3bda23f-6ec8-4977-9fad-e7cf99aa8bb3
# ╠═d5286246-3265-421e-b724-13b963f09887
# ╠═1420b89c-4bd5-4278-8a70-92304a9a5723
# ╠═e729ea16-387f-439d-843a-5383ad93703e
# ╟─5e7753e3-794a-4aa1-ae4a-a9cbedc9fc52
# ╠═fef94d49-27e4-4809-9804-270e9270f958
# ╠═5fff0d6b-5098-4500-a521-c0553ff5a4ef
# ╠═f0bd22ad-8e33-457e-b9b5-82805a3cc547
# ╠═0b396e36-3022-4c39-8e44-84ebd42d7b19
# ╠═9f0f1131-52cc-4894-a1e0-7999fbd5234a
# ╠═b2e2c302-eefd-41ec-97ca-aa302d62f2f1
# ╠═2cc30337-2c48-476f-b45f-c2bfb788a3db
# ╠═e7bb6538-e885-4ff9-8246-d888ced1da0f
# ╠═8330fc10-22c9-419b-887b-7d91ee580dc3
# ╠═cd469ae2-74d4-439e-ad27-9f2ea999a13c
