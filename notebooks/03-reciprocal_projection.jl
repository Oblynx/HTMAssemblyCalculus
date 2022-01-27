### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
begin
	import Pkg; Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using HierarchicalTemporalMemory
	using Random, Chain, Setfield, Statistics, Plots, PlutoUI
end

# ╔═╡ 8adb2e56-320c-420f-9938-3c96b7f06f22
"Include assembly operations library."
module lib
	include(joinpath(dirname(Base.active_project()), "src/assembly_operations.jl"))
end

# ╔═╡ 1bb0fcfc-2d7a-4634-9c93-263050c56a55
md"""
# Reciprocal projection

This notebook implements *reciprocal projection* between 2 HTM regions, building on the construction of the *projection* operation in [01-intro_projection.jl](./open?path=notebooks/01-intro_projection.jl).

Projection of assembly `x` from region A -> B is the operation where a stable assembly `y` is generated in B that fires consistently when `x` fires in A.
Reciprocal projection adds the requirement that *`x` fire in A when `y` fires in B.*

Instead of activating `x` externally, this will employ 3 Regions:

- input region `I`
- regions `A, B` with recurrent internal connections

The only difference between projection and reciprocal projection is the feedforward connection `B → A`

```
  I --> A --> B
        ↑     |
        +-----+
```

## Convergence

We will implement this simple algorithm with HTM regions and check the convergence of `x,y`.
If `xₜ and yₜ` both converge to a limit circle as $t→∞$ (practically 10-100 steps) of a small enough radius, we will consider the projection implemented.
"""

# ╔═╡ 6790a86a-3990-4c6a-878f-c18779f8d48d
md"""
Now all we need to do is to repeatedly stimulate `A` with `a`, propagating the stimulation of A->B and B->A according to the connection diagram.
In each timestep, `A` will produce a `xₜ` and `B` a `yₜ`:
"""

# ╔═╡ 2b6d959f-63bb-4d42-a34d-70d93f6a1636
md"""
The difference in activity is much less than `k` (neurons per minicolumn), but it is still evident.
This means that after convergence there are still many neurons firing in each minicolumn;
in fact, the average number of neurons firing in each active minicolumn should be close to the region's activity in the previous plot.
"""

# ╔═╡ dfa42c2b-6424-458a-a589-a1c5429ec4db
B= Region(lib.params_M.sp, lib.params_M.tm);

# ╔═╡ 043f8da8-f17a-4486-9c49-7d5ecc75457f
md"""
#### Regions

Let's define the HTM regions I, A, B with the same parameter values that were tuned for the projection operation.
"""

# ╔═╡ 354b27e3-75b4-45ce-9b45-88b8c110039f
I= Region(lib.params_input.sp, lib.params_input.tm);

# ╔═╡ 3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
a= @chain bitrand(lib.Nin) I(_).active;

# ╔═╡ 3a5e3aa0-e34e-457c-a405-9bc95cd88910
activity_per_active_minicolumn(y)= @chain y reshape(_,k,:) count(_, dims=1) filter(x->x>0, _) mean

# ╔═╡ 62d41be1-2970-48e1-b689-c2ecf1ab10ce
md"`y[T,experiments]` is the entire history of `yᵢ`"

# ╔═╡ 44343331-60aa-4985-b0c2-c7e9bc7885cd
HierarchicalTemporalMemory.gateCombine([a,a])|>typeof

# ╔═╡ bbaf0a64-2471-4106-a27c-9dd5a4f58c7c
A= Region(lib.params_M.sp, lib.params_M.tm)

# ╔═╡ 3df40867-674b-405a-a05c-0733a4caaa67
md"The 2 are almost identical, therefore this explanation of the first few steps in the graph seems to hold up."

# ╔═╡ a0971202-45c7-4947-a6f6-8ba1301d0316
reduce((a,b)-> HierarchicalTemporalMemory.bitVcat(a,b), [a,a])

# ╔═╡ 6a20715a-9dc9-461a-b6a2-f02522e277f0
md"We calculate the $Δyₜ=\texttt{reldistance}(yₜ,yₜ₊₁) ∀ t∈\{1..T-1\}$ and the experiment median. This metric will show whether `yₜ` converges."

# ╔═╡ bcba6521-1cfb-4ee7-a421-bd1def3770a4
T= 10; experiments= 1;

# ╔═╡ 7670ffbc-3621-4acb-b1e3-b56668470c24
y= Matrix(undef,T,experiments);

# ╔═╡ d08cce9d-dd9c-4530-82ae-bf1db8be3281
# Send SDRs from A -> B
with_terminal() do
	for e= 1:experiments
	  @info "experiment $e"
	  x= @chain bitrand(Nin) A(_).active;
	  reset!(B)
	  for t= 1:T
		t%10==0 && @info "t=$t"
		y[t,e]= step!(B,x).active
	  end
	end
end

# ╔═╡ e63a5860-587c-4aff-91be-a894b3443943
Δyₜₑ()= [reldistance(y[i,e], y[i+1,e]) for i=1:T-1, e=1:experiments];

# ╔═╡ 3229da2a-92f3-4064-bb80-cbd9f5523d7d
begin
	reset!(A)
	reset!(B)
	for t= 1:T
		xₜ= step!(A, [a,a])
	end
end

# ╔═╡ 56e0da0b-8858-459a-9aae-0eba4797cc06
md"""
## Understanding the convergence plot

As we see from the median (thick line), the experiments **converge after step $t=20$**.
This is double what assembly theory expects (10 steps), but it might be sensitive to the model's parameters.

So far, it looks like the "projection" operation has been implemented with HTM successfully.

#### Why is the distance exactly 0 in the beginning?

_Distal and proximal synapses_

This has to do with the biggest deviation of HTM from assembly calculus.
In assembly calculus, all synapses can trigger neurons to fire.
In HTM there are 2 kinds of synapses: _proximal_ and _distal_.
Only the proximal ones cause neurons to fire; distal synapses only cause neurons to win in competition among many that are proximally excited (keyword: *context*).

Neurons are grouped in "minicolumns", wherein they share the same proximal synapses, differing only in distal synapses.
Feedforward connectivity from A → B is through proximal synapses, while recurrent connectivity within B is through distal.
Every time `x` fires the same large set of neurons (minicolumns) is primed to fire in B.
In the beginning, before recurrent connections form through learning, every neuron in the primed minicolumns fires.
Gradually, recurrent connections form, which cause only a few neurons (often 1) in each minicolumn to be "excited" through distal connections.
But this distal excitement doesn't cause neurons to fire;
only, if the minicolumn happens to be activated in the following step, the distally excited neuron will win an in-column competition and be the only one to fire.

If this hypothesis is correct, then we expect to see a significant drop in the number of active neurons in B after the first few steps. Let's test this.
"""

# ╔═╡ 5e0ac810-1689-4bfb-88a8-85504cd821b6
md"""
### Repeated stimulation

To plot the convergence of $yᵢ$ we will now produce them repeatedly and store them.
Each experiment will run for `T` steps and we will run 15 experiments with different `x`.
"""

# ╔═╡ c866c71e-5dce-4af1-84bb-045815d1acb9
expmedian(y)= median(y, dims=2);

# ╔═╡ cf83549c-03b1-4f7d-be29-d3a29da3e8f7
activity_n(y)= count.(y)|> expmedian

# ╔═╡ 5ec611f4-d3ea-4c5d-9221-4dcf2da3e18c
abs.(expmedian(activity_per_active_minicolumn.(y)) .- activity_n(y)/activity_n(y)[1]*k)|>sum

# ╔═╡ 75560584-1631-4c93-8808-e269d345e1c8
Δȳₜ= Δyₜₑ()|> expmedian;

# ╔═╡ 85dfe3ed-cf70-4fde-a94a-f70c3fe4abf6
begin
	using Plots.PlotMeasures
	plot(activity_n(y)/activity_n(y)[1], linecolor=:blue, foreground_color_subplot=:blue, ylabel="Relative activity in B", rightmargin=14mm, legend=:none)
	plot!(twinx(), Δȳₜ, linecolor=:red, foreground_color_subplot=:red, ylabel="median Δȳ", legend=:none)
end

# ╔═╡ 456021aa-3dbf-4102-932f-43109905f9eb
begin
	plot(Δyₜₑ(), minorgrid=true, ylim= [0,.3], opacity=0.3, labels=nothing,
		title="Convergence of yₜ", ylabel="Relative distance per step Δyₜ",
		xlabel="t")
	plot!(Δȳₜ, linewidth=2, label="median Δȳ")
	hline!(zeros(T) .+ .05, linecolor=:black, linestyle=:dash, label="5% limit")
end

# ╔═╡ 91fa8319-1702-4cad-8ae8-0c1ed51a7bb8
md"""
This is the input activation `a`, which will stimulate the reciprocal projection:
"""

# ╔═╡ Cell order:
# ╠═0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
# ╠═8adb2e56-320c-420f-9938-3c96b7f06f22
# ╟─1bb0fcfc-2d7a-4634-9c93-263050c56a55
# ╠═6790a86a-3990-4c6a-878f-c18779f8d48d
# ╠═2b6d959f-63bb-4d42-a34d-70d93f6a1636
# ╠═7670ffbc-3621-4acb-b1e3-b56668470c24
# ╠═d08cce9d-dd9c-4530-82ae-bf1db8be3281
# ╠═cf83549c-03b1-4f7d-be29-d3a29da3e8f7
# ╠═5ec611f4-d3ea-4c5d-9221-4dcf2da3e18c
# ╠═dfa42c2b-6424-458a-a589-a1c5429ec4db
# ╠═85dfe3ed-cf70-4fde-a94a-f70c3fe4abf6
# ╟─043f8da8-f17a-4486-9c49-7d5ecc75457f
# ╠═354b27e3-75b4-45ce-9b45-88b8c110039f
# ╠═3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
# ╠═75560584-1631-4c93-8808-e269d345e1c8
# ╠═3a5e3aa0-e34e-457c-a405-9bc95cd88910
# ╠═62d41be1-2970-48e1-b689-c2ecf1ab10ce
# ╠═456021aa-3dbf-4102-932f-43109905f9eb
# ╠═44343331-60aa-4985-b0c2-c7e9bc7885cd
# ╠═e63a5860-587c-4aff-91be-a894b3443943
# ╠═bbaf0a64-2471-4106-a27c-9dd5a4f58c7c
# ╠═3df40867-674b-405a-a05c-0733a4caaa67
# ╠═3229da2a-92f3-4064-bb80-cbd9f5523d7d
# ╠═a0971202-45c7-4947-a6f6-8ba1301d0316
# ╟─6a20715a-9dc9-461a-b6a2-f02522e277f0
# ╠═bcba6521-1cfb-4ee7-a421-bd1def3770a4
# ╠═56e0da0b-8858-459a-9aae-0eba4797cc06
# ╠═5e0ac810-1689-4bfb-88a8-85504cd821b6
# ╠═c866c71e-5dce-4af1-84bb-045815d1acb9
# ╟─91fa8319-1702-4cad-8ae8-0c1ed51a7bb8
