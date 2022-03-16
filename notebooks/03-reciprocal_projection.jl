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
	using Images
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

We will activate `x` externally and employ 2 regions `A, B` projecting bidirectionally to each other.

The only difference between projection and reciprocal projection is the proximal ("feedforward") connection `B → A`

```
  in --> A --> B --> out
         ↑     |
         +-----+
```

## Convergence

We will implement this simple algorithm with HTM regions and check the convergence of `x,y`.
If `xₜ and yₜ` both converge to a limit circle as $t→∞$ (practically 10-100 steps) of a small enough radius, we will consider the projection implemented.

## Network to connect Regions

To implement this circuit diagram we're going to use an HTM `Network`, comprising the regions A,B and the adjacency matrix of their connections.
Let's create the regions and the network.
"""

# ╔═╡ 5420952d-4e11-469a-ab1e-67ed8fff6c4a
I= Region(lib.params_input.sp, lib.params_input.tm);

# ╔═╡ 1caa5db3-b761-43f9-bbb2-66c3b4779b56
A= Region(lib.params_M.sp, lib.params_M.tm);

# ╔═╡ e7fb5232-502d-4c94-8250-4fff4b7b4176
B= Region(lib.params_M.sp, lib.params_M.tm);

# ╔═╡ becb08ed-fbec-4047-b8d3-a29cae874114
md"""
The network takes a vector with the regions A,B and an adjacency matrix with their connections.
The circuit diagram
```
  in --> A --> B --> out
         ↑     |
         +-----+
```

would suggest the adjacency matrix

$(Gray.([
	0 1
	1 0
]))

However, the network also needs to connect 2 extra "virtual" regions, the input and output. This works by appending 2 extra elements to the end of the adjacency matrix, connecting the 3rd row (input) to region 1 (A) and connecting region 2 (B) to the output:

$(Gray.([
	0 1 0 0
    1 0 0 1
	1 0 0 0
	0 0 0 0
]))

This is the adjacency matrix that we'll use for the network.
"""

# ╔═╡ 2fedb87d-85ab-4c11-8e6c-2494a3311b38
reciprocal_net= Network([A,B], connect_forward= [
	0 1 0 0
	1 0 0 1
	1 0 0 0
	0 0 0 0
]);

# ╔═╡ 6790a86a-3990-4c6a-878f-c18779f8d48d
md"""
Now all we need to do is to repeatedly stimulate `A` with `a`, propagating the stimulation of A->B and B->A according to the connection diagram.
In each timestep, `A` will produce a `xₜ` and `B` a `yₜ`:
"""

# ╔═╡ 043f8da8-f17a-4486-9c49-7d5ecc75457f
md"""
## Stimulating the network

Stimulating the network with the familiar `step!` function means to present the stimulus at the input node and let each region process the inputs presently available to it at each edge of the network.

Presenting `a` at the input node with `reciprocal_net(a)` is going to be equivalent to presenting `a` to region A, or `A(a)`. However, in the case of the network we will also calculate the activation of `B` with the input that is available to it, which is initially 0.

With each step, the signal is propagating through the network by 1 edge. If we present again `a` at step 2, the external input of `A` will stay the same, but the `B` will now be activated by `A(a)`.
Taking the reciprocal connection into account, the recurrent equations of the network for time $t$ are as follows.
With $A$ we symbolize the activation function of a region and $Aₜ$ is its output at time $t$.
$a|b$ means bitwise OR. 

``Aₜ = A(a | Bₜ₋₁)``

``Bₜ = B(Aₜ₋₁)``

The same propagation rules apply for the `step!` function, which simply calls `step!` on each region.

#### Experiment
"""

# ╔═╡ 3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
a()= @chain bitrand(lib.Nin) I(_).active;

# ╔═╡ bcba6521-1cfb-4ee7-a421-bd1def3770a4
T= 30; experiments= 5; net_T= 4T

# ╔═╡ 57f45bc3-4978-4abc-8763-0bcd9ce62542
md"""
Repeatedly stimulate the network with `a`. If `T` is the latency of the projection operation, an estimated time for convergence should be around `2T` (in fact, this is seen when simulating a feedforward network with 2 regions).
We will capture the activation of each region over time and plot their convergence.
We will also capture the interconnectivity and other metrics used in the Projection notebook.

Let's collect the data:
"""

# ╔═╡ 274fe051-0317-4651-b6eb-7c7824b7c5b1
begin
	# activation[t][region] and connectivity metrics
	stepRepeatedly!(net,a)= [ begin
		step!(net,a)
		(
			activity= net.region_α₋[:],
			interconnected= lib.interconnectionMeasure.(net.region_α₋, net.regions),
			density= lib.interconnectionDensity.(net.region_α₋, net.regions),
		)
	end for t= 1:net_T ]
	experiment(id,a)= stepRepeatedly!(deepcopy(reciprocal_net), a)
	reset!(reciprocal_net)
	α= @chain begin
	  @sync [Threads.@spawn experiment(e,a()) for e= 1:experiments]
	  fetch.(_)
	end
end

# ╔═╡ d4fa5204-79ab-4dd0-88be-92b65a546580
md"""
The experiment results `α` have 3 indices and a selector: `α[experiment][time]{activity,robustness}[region]`.
Let's split them into 2 variables for each region, `A -> x`, `B -> y` with indices `x[experiment][time]`.
"""

# ╔═╡ b4b282ee-83d0-41ce-8f7c-eca27e13b78a
parseMeasurements(x,field,region)= [ map(e-> map(t-> getproperty(t,field)[region], e),x)[e][t] for t=1:net_T, e=1:experiments];

# ╔═╡ 287e2f5d-22a6-40eb-a534-2e0c1de744b2
begin
	x= parseMeasurements(α,:activity,1);
	x_interconnected= parseMeasurements(α,:interconnected,1);
	x_density= parseMeasurements(α,:density,1);
	y= parseMeasurements(α,:activity,2);
	y_interconnected= parseMeasurements(α,:interconnected,2);
	nothing
end

# ╔═╡ addb6923-725a-451e-a1fb-af35c3cd6aee
md"""
The convergence curve:
"""

# ╔═╡ 396e1062-1acc-4aa0-ac29-fa245b0f2f46
Δy= [lib.reldistance(y[i,e], y[i+1,e]) for i=1:net_T-1, e=1:experiments];

# ╔═╡ 5e4444cd-41bc-4720-a18c-657c8b7ee3ee
Δx= [lib.reldistance(x[i,e], x[i+1,e]) for i=1:net_T-1, e=1:experiments];

# ╔═╡ ce3ded6c-889e-4db1-b40c-0548be4397fd
expmedian(y)= median(y, dims=2);

# ╔═╡ befd54da-edf0-4fd2-9cbd-752fed1fd60c
begin
	plot(Δy, minorgrid=true, ylim= [0,1], opacity=0.3, labels=nothing,
		title="Convergence of yₜ", ylabel="Relative distance per step Δyₜ",
		xlabel="t")
	plot!(Δy|>expmedian, linewidth=2, label="median Δȳ")
	hline!(zeros(T) .+ .05, linecolor=:black, linestyle=:dash, label="5% limit")
end

# ╔═╡ c3f555dc-e80b-4d2c-9969-fd617219aee3
normalize(v; dims=1)= (v .- minimum(v,dims=dims)) ./ maximum(v,dims=dims)

# ╔═╡ 2a119b92-d2e2-4799-ae32-a66f00d4edd4
begin
	plot(hcat(x_interconnected|> expmedian|> normalize,
		y_interconnected|> expmedian|> normalize),
		minorgrid=true, title="Assembly interconnection density training curve",
		ylabel="interconnection density",
		xlabel="t", label=:none
	)
	vline!([30,30],linecolor=:black, opacity=0.3, linestyle=:dash, label=:none)
end

# ╔═╡ 1a1c1ad8-5266-4f70-87f5-adb88828b1e6
md"""
## Pretraining of `x`

Create first assembly `x` in `A`, then train network.
"""

# ╔═╡ 49baa0cc-cb74-4f0d-a966-b153013a649f
begin
	reset!(reciprocal_net)
	experimentPre(id)= begin
		net= deepcopy(reciprocal_net)
		a_pre= a()
		lib.project!(net.regions[1],a_pre)
		stepRepeatedly!(net, a_pre)
	end
	α_pretrain= @chain begin
	  @sync [Threads.@spawn experimentPre(e) for e= 1:experiments]
	  fetch.(_)
	end
end

# ╔═╡ ff5993c2-8705-4e62-afe7-5aab0549fc18
begin
	x_pre= parseMeasurements(α_pretrain,:activity,1);
	x_pre_interconnected= parseMeasurements(α_pretrain,:interconnected,1);
	x_pre_density= parseMeasurements(α_pretrain,:density,1);
	y_pre= parseMeasurements(α_pretrain,:activity,2);
	y_pre_interconnected= parseMeasurements(α_pretrain,:interconnected,2);
	nothing
end

# ╔═╡ 867d0540-62a1-4675-9ef5-8e712b1a82c4
Δx_pre= [lib.reldistance(x_pre[i,e], x_pre[i+1,e]) for i=1:net_T-1, e=1:experiments];

# ╔═╡ 9e06b199-b05b-4c0b-b075-a98e526704c2
Δy_pre= [lib.reldistance(y_pre[i,e], y_pre[i+1,e]) for i=1:net_T-1, e=1:experiments];

# ╔═╡ 9e7f3181-6de9-4cf9-b536-e82f8d02db9f
begin
	plot(Δx_pre, minorgrid=true, ylim= [0,1], opacity=0.3, labels=nothing,
		title="Convergence of yₜ", ylabel="Relative distance per step Δyₜ",
		xlabel="t")
	plot!(Δx_pre|>expmedian, linewidth=2, label="median Δȳ")
	hline!(zeros(T) .+ .05, linecolor=:black, linestyle=:dash, label="5% limit")
end

# ╔═╡ 998f3cb3-1098-4262-8216-d5f48b72d6bd
begin
	plot(hcat(x_pre_interconnected|> expmedian|> normalize,
		y_pre_interconnected|> expmedian|> normalize),
		minorgrid=true, title="Assembly interconnection density training curve",
		ylabel="interconnection density",
		xlabel="t", label=:none
	)
	vline!([30,30],linecolor=:black, opacity=0.3, linestyle=:dash, label=:none)
end

# ╔═╡ Cell order:
# ╠═0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
# ╠═8adb2e56-320c-420f-9938-3c96b7f06f22
# ╟─1bb0fcfc-2d7a-4634-9c93-263050c56a55
# ╠═5420952d-4e11-469a-ab1e-67ed8fff6c4a
# ╠═1caa5db3-b761-43f9-bbb2-66c3b4779b56
# ╠═e7fb5232-502d-4c94-8250-4fff4b7b4176
# ╟─becb08ed-fbec-4047-b8d3-a29cae874114
# ╠═2fedb87d-85ab-4c11-8e6c-2494a3311b38
# ╟─6790a86a-3990-4c6a-878f-c18779f8d48d
# ╟─043f8da8-f17a-4486-9c49-7d5ecc75457f
# ╠═3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
# ╠═bcba6521-1cfb-4ee7-a421-bd1def3770a4
# ╟─57f45bc3-4978-4abc-8763-0bcd9ce62542
# ╠═274fe051-0317-4651-b6eb-7c7824b7c5b1
# ╠═d4fa5204-79ab-4dd0-88be-92b65a546580
# ╠═b4b282ee-83d0-41ce-8f7c-eca27e13b78a
# ╠═287e2f5d-22a6-40eb-a534-2e0c1de744b2
# ╟─addb6923-725a-451e-a1fb-af35c3cd6aee
# ╠═396e1062-1acc-4aa0-ac29-fa245b0f2f46
# ╠═5e4444cd-41bc-4720-a18c-657c8b7ee3ee
# ╠═ce3ded6c-889e-4db1-b40c-0548be4397fd
# ╠═befd54da-edf0-4fd2-9cbd-752fed1fd60c
# ╠═c3f555dc-e80b-4d2c-9969-fd617219aee3
# ╠═2a119b92-d2e2-4799-ae32-a66f00d4edd4
# ╠═1a1c1ad8-5266-4f70-87f5-adb88828b1e6
# ╠═49baa0cc-cb74-4f0d-a966-b153013a649f
# ╠═ff5993c2-8705-4e62-afe7-5aab0549fc18
# ╠═867d0540-62a1-4675-9ef5-8e712b1a82c4
# ╠═9e06b199-b05b-4c0b-b075-a98e526704c2
# ╠═9e7f3181-6de9-4cf9-b536-e82f8d02db9f
# ╠═998f3cb3-1098-4262-8216-d5f48b72d6bd
