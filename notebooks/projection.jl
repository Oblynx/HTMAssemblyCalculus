### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
begin
	import Pkg; Pkg.activate(Base.current_project())
	using Revise
	using HierarchicalTemporalMemory
	using Random, Chain, Setfield, Statistics, Plots, PlutoUI
	#using PlotlyJS;	plotlyjs()
end

# ╔═╡ 1bb0fcfc-2d7a-4634-9c93-263050c56a55
md"""
# Projection

This notebook implements the simplest operation of assembly calculus, *projection*, with HTM.

Assembly calculus and HTM share the concept of "Region" as a modular circuit of neurons.
HTM imposes specific structure in the Region.
The question is: _is the HTM structure behaving according to assembly calculus_?

Projection of assembly `x` from region A -> B is the operation where a stable assembly `y` is generated in B that fires consistently when `x` fires in A.
In assembly calculus this is shown to happen simply by repeated firing of `x` for a number of time steps:

```algo
repeat T:    # T: constant ~10
  fire x
```

Alternatively (to not require a constant T):

```algo
repeat until yₜ converges:
  fire x
```

Note that we don't necessarily expect point stability, but a limit circle of low enough radius. More on this with the convergence measure.

## Experiment

### The Plan

We will implement this simple algorithm with HTM regions and check the convergence of `y`.
If `yₜ` converges to a limit circle as $t→∞$ (practically 10-100 steps) of a small enough radius, we will consider the projection implemented.

But how small is small enough?
The algorithm parameters might have a big impact on the convergence radius, or might even prohibit convergence.
Intuitively in order to set an expectation, let the limit radius $|yₜ₁ - yₜ₂| < 5\%|y|$.

Let's define the relative distance behind this limit radius.
Since $yᵢ$ are sparse binary vectors:
- the core "similarity" kernel will be bitwise AND: $yᵢ .\& yⱼ$
- the norm $|yᵢ|$ will be the count of `true`

The relative distance $\frac{|yᵢ - yⱼ|}{|y|}$ can be:
"""

# ╔═╡ 70471bdc-660e-442e-b92a-f486abd120a7
reldistance(yᵢ,yⱼ)= 1 - count(yᵢ .& yⱼ) / min(count(yᵢ), count(yⱼ))

# ╔═╡ 7620c202-c6cf-44db-9cda-13642e28b45a
md"""
which is bounded by $0 ≤ \texttt{reldistance}(yᵢ, yⱼ) ≤ 1$
"""

# ╔═╡ 043f8da8-f17a-4486-9c49-7d5ecc75457f
md"""
### Collect the ingredients

Let's define the HTM regions A, B with the same properties (adapting the input size to match):
"""

# ╔═╡ bbaf0a64-2471-4106-a27c-9dd5a4f58c7c
begin
	Nin= 1e3|> Int
	Nn= 4e5
	k= 10
	_Nc()= floor(Int,Nn/k)
	sparsity()= 1/sqrt(_Nc())
	params_A= (
		sp= SPParams(szᵢₙ=Nin, szₛₚ=_Nc(), s= sparsity(), prob_synapse=1e-3, enable_local_inhibit=false),
		tm= TMParams(Nc=_Nc(), k=k, p⁺_01= .11, θ_stimulus_learn=15, θ_stimulus_activate= 22)
	)
	params_B= @set params_A.sp.szᵢₙ= params_A.tm.Nₙ
end

# ╔═╡ 334b271d-247c-465c-ae82-f91007a6d9d0
md"The regions:"

# ╔═╡ 4833f12f-3eac-407c-a9ce-0d1aea216077
begin
	A= Region(params_A.sp, params_A.tm)
	B= Region(params_B.sp, params_B.tm)
	A,B
end

# ╔═╡ 91fa8319-1702-4cad-8ae8-0c1ed51a7bb8
md"""
This is the assembly `x`, which first needs a random input.
"""

# ╔═╡ 3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
x= @chain bitrand(Nin) A(_).active;

# ╔═╡ 6790a86a-3990-4c6a-878f-c18779f8d48d
md"""
Now all we need to do is to repeatedly stimulate `B` with `x` and let it learn.
Each time `B` is stimulated, it will produce a `yᵢ`:
"""

# ╔═╡ 5458c6ae-ace2-4928-8a9b-ef919a1e97cb
y₁= B(x).active;

# ╔═╡ 6b56422c-5f9b-41d8-a080-ef5ca3d7db7b
md"""
But in order to let region `B` learn from `x`, we need to call the `step!` function, which will return the same output as `yᵢ` the first time it's called, while at the same time adapting synapses:
"""

# ╔═╡ 3229da2a-92f3-4064-bb80-cbd9f5523d7d
begin
	reset!(B)
	(step!(B,x).active .== y₁) |> all
end

# ╔═╡ 5e0ac810-1689-4bfb-88a8-85504cd821b6
md"""
### Repeated stimulation

To plot the convergence of $yᵢ$ we will now produce them repeatedly and store them.
Each experiment will run for `T` steps and we will run 10 experiments with different `x`.
`y` is the entire history of `yᵢ`
"""

# ╔═╡ 1015f629-9818-4dbb-b574-97c552e96164
begin
	T= 40; experiments=15;
	y= Matrix(undef,T,experiments);
end

# ╔═╡ d08cce9d-dd9c-4530-82ae-bf1db8be3281
# Send SDRs from A -> B
for e= 1:experiments
  x= @chain bitrand(Nin) A(_).active;
  reset!(B)
  for t= 1:T
	#t%10==0 && @info "t=$t"
	y[t,e]= step!(B,x).active
  end
end

# ╔═╡ e63a5860-587c-4aff-91be-a894b3443943
similarity= [reldistance(y[i,e], y[i+1,e]) for i=1:T-1, e=1:experiments];

# ╔═╡ 456021aa-3dbf-4102-932f-43109905f9eb
begin
	plot(similarity, minorgrid=true, ylim= [0,.3], opacity=0.3, legend=nothing,
		title="Convergence of yᵢ", ylabel="Relative distance",
		xlabel="Steps")
	plot!(median(similarity,dims=2), linewidth=2, label="ȳ")
	hline!(zeros(size(similarity,1)) .+ .05, linecolor=:black, linestyle=:dash)
end

# ╔═╡ Cell order:
# ╟─0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
# ╟─1bb0fcfc-2d7a-4634-9c93-263050c56a55
# ╠═70471bdc-660e-442e-b92a-f486abd120a7
# ╟─7620c202-c6cf-44db-9cda-13642e28b45a
# ╟─043f8da8-f17a-4486-9c49-7d5ecc75457f
# ╠═bbaf0a64-2471-4106-a27c-9dd5a4f58c7c
# ╟─334b271d-247c-465c-ae82-f91007a6d9d0
# ╠═4833f12f-3eac-407c-a9ce-0d1aea216077
# ╟─91fa8319-1702-4cad-8ae8-0c1ed51a7bb8
# ╠═3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
# ╟─6790a86a-3990-4c6a-878f-c18779f8d48d
# ╠═5458c6ae-ace2-4928-8a9b-ef919a1e97cb
# ╟─6b56422c-5f9b-41d8-a080-ef5ca3d7db7b
# ╠═3229da2a-92f3-4064-bb80-cbd9f5523d7d
# ╟─5e0ac810-1689-4bfb-88a8-85504cd821b6
# ╠═1015f629-9818-4dbb-b574-97c552e96164
# ╠═d08cce9d-dd9c-4530-82ae-bf1db8be3281
# ╠═e63a5860-587c-4aff-91be-a894b3443943
# ╠═456021aa-3dbf-4102-932f-43109905f9eb
