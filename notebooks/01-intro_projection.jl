### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
begin
	import Pkg; Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using HierarchicalTemporalMemory
	using Random, Chain, Setfield, Statistics, Plots, PlutoUI
	using Plots.PlotMeasures
end

# ╔═╡ 1bb0fcfc-2d7a-4634-9c93-263050c56a55
md"""
# What is Assembly calculus?

Neuron Assembly Calculus is a computational framework for cognitive function. It describes a dynamical system of brain areas with random connections between excitatory neurons, with Hebbian plasticity, where only the top-k most activated neurons fire. From these properties **neuron assemblies** emerge: clusters of highly interconnected neurons in the same area, which can be created through programmatic operations. The authors probabilistically prove the convergence of these operations based on the model's structure and connectome (ref: [Neuron assembly calculus](https://www.pnas.org/content/117/25/14464), [A biologically plausible parser](https://direct.mit.edu/tacl/article/doi/10.1162/tacl_a_00432/108608/A-Biologically-Plausible-Parser)).

## Hierarchical Temporal Memory

Hierarchical Temporal Memory (HTM) is a biologically constrained model of brain computation. It too is a dynamical system of brain areas with excitatory neurons, Hebbian plasticity and top-k activation. However, the dynamics are more strongly constrained. For example, neurons are arranged in local circuits with local competition, and have multiple dendrites, where they receive less-stimulating input ([Why neurons have thousands of synapses](https://www.frontiersin.org/articles/10.3389/fncir.2016.00023/full)). This model of a brain area's dynamics underpins the [Thousand Brains theory](https://link.springer.com/article/10.1007/s42452-021-04715-0) that attempts to explain human neocortical computation.

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

Note that we don't necessarily expect point stability, but that no new neurons are being recruited to represent the output.
Point stability would correspond to getting exactly the same output yₜ ∀ t>tₛ.
Instead, convergence of yₜ relies on the set of all neurons that were activated until now (ie. belonging in any yₜ ∀ t<tₛ).

### Convergence metric

Let Sₜ be the set of neurons that have been activated in the region between time 0 and t:
``Sₜ = \sum_0^t y₀ + y₁ + ... + yₜ``,
where ``+`` is boolean OR.
Convergence is achieved at t₀ when ``Sₜ = S\\_{t-1} ∀ t>t₀`` .

In practice, we want a limited convergence metric that only looks back a set amount of steps `τ`, not to the beginning of time.
This `S(yₜ,τ)` is a moving average filter on yₜ.
"""

# ╔═╡ 2e40a09d-3220-4424-851d-8470dfa235d1
S(y,τ)= [reduce((a,b)->a.|b, @view y[max(1,t-τ):t]) for t=1:length(y)]

# ╔═╡ e109962a-dfaa-438d-aa86-3e074e57f51d
md"""
## Experiment

### The Plan

We will implement this simple algorithm with HTM regions and check the convergence of `y` using the convergence metric `S`.
If `Sₜ` converges to a limit circle as $t→∞$ (practically 10-100 steps) of a small enough radius, we will consider the projection implemented.

But how small is small enough?
The algorithm parameters might have a big impact on the convergence radius, or might even prohibit convergence.
Intuitively in order to set an expectation, let the limit radius $|Sₜ₁ - Sₜ₂| < 5\%|S|$.

Let's define the relative distance behind this limit radius.
Since $yᵢ$ are sparse binary vectors:
- the core "similarity" kernel will be bitwise AND: $yᵢ \& yⱼ$
- the norm $|yᵢ|$ will be the count of `true`
- the similarity will be normalized by the harmonic mean of $|yᵢ|,|yⱼ|$

The relative distance $\frac{|yᵢ - yⱼ|}{|y|}$ can be:
"""

# ╔═╡ 70471bdc-660e-442e-b92a-f486abd120a7
reldistance(yᵢ,yⱼ)= 1 - count(yᵢ .& yⱼ) * mean(inv,[count(yᵢ),count(yⱼ)])

# ╔═╡ 7620c202-c6cf-44db-9cda-13642e28b45a
md"""
which is bounded by $0 ≤ \texttt{reldistance}(yᵢ, yⱼ) ≤ 1$
"""

# ╔═╡ 043f8da8-f17a-4486-9c49-7d5ecc75457f
md"""
### Collect the ingredients

Let's define the HTM regions A, B with the same properties (adapting the input size to match).
Minicolumns will be explained later.
"""

# ╔═╡ bbaf0a64-2471-4106-a27c-9dd5a4f58c7c
begin
	Nin= 1e3|> Int        		# input size
	Nn= 20e4             	    # number of neurons in each area
	k= 15                 		# neurons per minicolumn
	thresholds= ( 				# calibrate how many dendritic inputs needed to fire
		tm_learn= 14,
		tm_activate= 19,
		tm_dendriteSynapses= 60,
	)
	learnrate= (
		dist_p⁺= .058,
		dist_p⁻= .015,
		dist_LTD_p⁻= .0001,
		prox_p⁺= .10,
		prox_p⁻= .04,
	)
end

# ╔═╡ 9c7c7c2b-b020-4a83-8917-909470331be7
begin
	_Nc()= floor(Int,Nn/k) 		# number of minicolumns
	sparsity()= 1/sqrt(_Nc());
	params_A= (
		sp= SPParams(szᵢₙ=Nin, szₛₚ=_Nc(), s= sparsity(), prob_synapse=6e-3,
			p⁺_01= learnrate.prox_p⁺, p⁻_01= learnrate.prox_p⁻,
			enable_local_inhibit=false),
		tm= TMParams(Nc=_Nc(), k=k,
			p⁺_01= learnrate.dist_p⁺, p⁻_01= learnrate.dist_p⁻, LTD_p⁻_01= learnrate.dist_LTD_p⁻,
			θ_stimulus_learn=thresholds.tm_learn,
			θ_stimulus_activate=thresholds.tm_activate,
			synapseSampleSize=thresholds.tm_dendriteSynapses)
	)
	# parameters for "merge" or "work" region, where assemblies will be formed
	params_B= @set params_A.sp.szᵢₙ= params_A.tm.Nₙ
end

# ╔═╡ 334b271d-247c-465c-ae82-f91007a6d9d0
md"The regions:"

# ╔═╡ 4833f12f-3eac-407c-a9ce-0d1aea216077
A= Region(params_A.sp, params_A.tm);

# ╔═╡ 45d8e50b-8baa-4da2-8cd0-5c81d634c450
B= Region(params_B.sp, params_B.tm);

# ╔═╡ 91fa8319-1702-4cad-8ae8-0c1ed51a7bb8
md"""
This is the assembly `x`, which first needs a random input.
"""

# ╔═╡ 3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
x= @chain bitrand(Nin) A(_).active;

# ╔═╡ 6790a86a-3990-4c6a-878f-c18779f8d48d
md"""
### The Projection operation

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
	md"`step!(B,x).active == y₁`: $((step!(B,x).active .== y₁) |> all)"
end

# ╔═╡ 71f21302-04db-41aa-a259-f3be2b7bc271
md"We will now explore the result of repeated stimulation"

# ╔═╡ 5e0ac810-1689-4bfb-88a8-85504cd821b6
md"""
### Repeated stimulation

To plot the convergence of $yᵢ$ we will now produce them repeatedly and store them.
Each experiment will run for `T` steps and we will run many experiments with different `x` each.
"""

# ╔═╡ 1015f629-9818-4dbb-b574-97c552e96164
T= 60; experiments=6;

# ╔═╡ 62d41be1-2970-48e1-b689-c2ecf1ab10ce
md"""
The results will be:

- `y[T,experiments]` is the entire history of `yᵢ`
- `yState[experiments]` is the end state of the region B for each experiment
"""

# ╔═╡ d08cce9d-dd9c-4530-82ae-bf1db8be3281
# Send SDRs from A -> B
begin
	randomX()= @chain bitrand(Nin) A(_).active;
	stepRepeatedly(B,x)= [step!(B,x).active for t= 1:T], B
	experiment(id)= stepRepeatedly(
		deepcopy(B),
		randomX()
	)
	reshapeMatrix(x,d1,d2)= @chain x begin
	  Iterators.flatten
      collect
      reshape(_, d1,d2)
	end
	splitReshapeExp(experimentVec)= (
		( @chain experimentVec first.(_) reshapeMatrix(_,T, experiments) ),
		( @chain experimentVec Base.rest.(_,2) first.(_) ),
	)
	reset!(B)
	y, yState= @chain begin
	  @sync [Threads.@spawn experiment(e) for e= 1:experiments]
	  fetch.(_)
	  splitReshapeExp
	end
	nothing
end

# ╔═╡ 6a20715a-9dc9-461a-b6a2-f02522e277f0
md"We calculate the $ΔSₜ=\texttt{reldistance}(Sₜ,Sₜ₊₁) ∀ t∈\{1..T-1\}$ and the experiment median. This metric will show whether `yₜ` converges."

# ╔═╡ 18e6d1c5-cd7a-4c44-a186-3ec2404f4ffc
Δrel(x)= @views [0; map(reldistance, x[1:end-1], x[2:end])]

# ╔═╡ d59c79f9-c570-4f71-a379-6808bc189e02
convergence(yₜ; τ=5)= Δrel(S(yₜ,τ))

# ╔═╡ e63a5860-587c-4aff-91be-a894b3443943
ΔSₜₑ()= @chain begin
	[convergence(y[:,e],τ=5) for e=1:experiments]
	Iterators.flatten
	collect
	reshape(_, :,experiments)
end

# ╔═╡ c866c71e-5dce-4af1-84bb-045815d1acb9
expmedian(y)= median(y, dims=2);

# ╔═╡ 75560584-1631-4c93-8808-e269d345e1c8
ΔŜₜ= ΔSₜₑ()|> expmedian;

# ╔═╡ 456021aa-3dbf-4102-932f-43109905f9eb
begin
	plot(ΔSₜₑ(), minorgrid=true, ylim= [0,.5], opacity=0.3, labels=nothing,
		title="Convergence of yₜ", ylabel="Relative distance ΔSₜ",
		xlabel="t")
	plot!(ΔŜₜ, linewidth=2, label="median ΔŜ")
	hline!(zeros(T) .+ .05, linecolor=:black, linestyle=:dash, label="5% limit")
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
Only the proximal ones cause neurons to fire; distal synapses only cause neurons to win in competition among many that are proximally excited in the same local circuit.

Neurons are grouped in "minicolumns", wherein they share the same proximal synapses, differing only in distal synapses.
Feedforward connectivity from A → B is through proximal synapses, while recurrent connectivity within B is through distal.
Every time `x` fires the same large set of neurons (minicolumns) is primed to fire in B.
In the beginning, before recurrent connections form through learning, every neuron in the primed minicolumns fires. In this case the minicolumn is **bursting**.
Gradually, recurrent connections form, which cause only a few neurons (often 1) in each minicolumn to be "excited" through distal connections.
But this distal excitement doesn't cause neurons to fire;
only, if the minicolumn happens to be activated in the following step, the distally excited neuron will win an in-column competition and be the only one to fire.

If this hypothesis is correct, then we expect to see a significant drop in the number of active neurons in B after the first few steps. Let's test this.
"""

# ╔═╡ cf83549c-03b1-4f7d-be29-d3a29da3e8f7
activity_n(y)= count.(y)|> expmedian

# ╔═╡ 85dfe3ed-cf70-4fde-a94a-f70c3fe4abf6
begin
	plot(activity_n(y)/activity_n(y)[1], ylim=[0,1],
		linecolor=:blue, xlabel="t", ylabel="Fraction of active neurons",
		rightmargin=14mm, legend=:none, title="Less surprise → fewer active neurons",
		minorgrid= true,
	)
	plot!(twinx(), ΔŜₜ, linecolor=:red,
		foreground_color_border=:red, foreground_color_axis=:red,
		ylabel="Assembly convergence", legend=:none
	)
end

# ╔═╡ 71a2a9e5-f5cd-4c61-bfc7-4ee950897e16
md"""
At the final moment of the simulation we can observe:
- Ideal fraction of active neurons: **$(1/k*100)%** (ie. ``1/k`` for perfect prediction)
- Actual fraction of active neurons: **$(round(activity_n(y)[end]/activity_n(y)[1]*100, sigdigits=2))%**
"""

# ╔═╡ 2b6d959f-63bb-4d42-a34d-70d93f6a1636
md"""
The difference in activity is much less than `k` (neurons per minicolumn), but it is still evident.
This means that after convergence there are still many neurons firing in each minicolumn;
in fact, the average number of neurons firing in each active minicolumn should be close to the region's activity in the previous plot.
"""

# ╔═╡ 3a5e3aa0-e34e-457c-a405-9bc95cd88910
activity_per_active_minicolumn(y)= @chain y reshape(_,k,:) count(_, dims=1) filter(x->x>0, _) mean

# ╔═╡ 5ec611f4-d3ea-4c5d-9221-4dcf2da3e18c
abs.(expmedian(activity_per_active_minicolumn.(y)) .- activity_n(y)/activity_n(y)[1]*k)|>sum

# ╔═╡ 3df40867-674b-405a-a05c-0733a4caaa67
md"The 2 are almost identical, therefore this explanation of the first few steps in the graph seems to hold up."

# ╔═╡ c3896446-3a24-413a-afab-6a67b07092f5
md"""
## Assemblies have dense interconnections

Assemblies are expected to have denser interconnections than random collections of neurons. To measure this property, we define the interconnection measure to be the number of synapses with pre- and post-synaptic neurons in the same population:
"""

# ╔═╡ dcbe4617-7314-42c1-87bf-57ca82554a20
interconnectionMeasure(activation, region)= activation' * distalSynapses(region) * activation

# ╔═╡ 059f3e13-56f6-4dc6-9570-19c74996d1ef
md"""
Using the interconnection measure we can now determine the average interconnectivity of the assemblies across experiments and compare with the average interconnectivity of 30 random neuronal activations that aren't assemblies.
"""

# ╔═╡ 8e9003fa-c822-4508-9862-58a816d9242d
assemblyInterconnection= map(interconnectionMeasure, y[end,:], yState)|> mean

# ╔═╡ d4cfbfe7-ba6f-4936-8854-8d49277657d2
randomInterconection= map(yState) do (r)
	@chain begin
		[bitrand(Nₙ(r)) for i= 1:30]       # 30 random activations
		interconnectionMeasure.(_,Ref(r))  # interconnection measure for each
		mean                               # mean interconnection among the 30
	end
end |> mean                                # mean across the experiments

# ╔═╡ dfff3823-d750-40e1-9ddf-769cb2c7e575
md"The average interconnectivity of assemblies is **×$(round(assemblyInterconnection / randomInterconection, digits=1))** times higher than random. It probably increases as the training process continues. Let's confirm this with a training graph:"

# ╔═╡ ed61ded0-d665-46a9-bed8-cd83cb0a545a
begin
	train_interconnection!(R,x)= begin
	    y= step!(R,x).active
		interconnectionMeasure(y,R)
	end
	Random.seed!(0)
	reset!(B)
	trainingcurve_interconnection= [train_interconnection!(B,x) for t= 1:1.5T]
	md"This cell produces the training curve."
end

# ╔═╡ 5bd5049c-de0d-4838-9c82-9cef604e650e
begin
	plot(trainingcurve_interconnection,
		minorgrid=true, title="Assembly interconnection training curve",
		ylabel="assembly interconnection",
		xlabel="t", label=:none
	)
	vline!([30,30],linecolor=:black, opacity=0.3, linestyle=:dash, label=:none)
end

# ╔═╡ b886095d-d449-4334-99cf-44b800bb4fc4
md"In the graph we can see that the interconnection measure reaches the limit cycle close to the cutoff $T=30$ that we chose to consider assemblies stable"

# ╔═╡ 8fdacada-55df-4abd-9501-7405e608529b
md"""
#### Is `x` an assembly?
### Note on defining assemblies and on input regions

In this notebook, we activate region `A` arbitrarily, through an external mechanism.
Is the resulting activation `x` an assembly?

The interconnection measure of `x` is $(interconnectionMeasure(x,A)).
That's because region A didn't go through any learning steps and regions begin with no recurrent connections; they only build them as needed by the learning rules.
`x` is a neuronal activation, but strictly speaking not an assembly.

However, this is not actually an important question.
An HTM region expects external stimulation. In this case, the external stimulation for B comes from A.
Correspondingly, we could have used another region to produce an external stimulation for A.
As long as we can assume:

1. some external stimulation
1. convergence of projection

we can use "input regions" that produce neuronal stimulations as input to others, without having to apply the learning rules themselves.
"""

# ╔═╡ 6a5298f1-dc8a-4422-b034-70df607c0456
md"""
## Assemblies as clusters

A cluster is a set of nodes with higher interconnection than their connection to outside elements.
We can define a new measure to probe the "clustering" of assemblies by comparing the interconnections to the external connections, that will be called "interconnectionDensity":
"""

# ╔═╡ 16fa7cce-a1db-4ffd-b769-3f0109305590
"`interconnectionDensity(activation, region)` is the ratio of connected synapses whose pre- and post-synaptic neurons are inside the activation set versus those whose pre- synaptic neurons are outside it."
interconnectionDensity(activation, region)= activation' * distalSynapses(region) * activation / ((.!activation)' * distalSynapses(region) * activation)

# ╔═╡ e9a82f4d-1866-45a2-a7e0-a5e08ce5e0fb
begin
	train_interconnectionDensity!(R,x)= begin
	    y= step!(R,x).active
		interconnectionDensity(y,R)
	end
	Random.seed!(0)
	reset!(B)
	trainingcurveDensity= [train_interconnectionDensity!(B,x) for t= 1:1.5T]
	md"This cell produces the training curve."
end

# ╔═╡ 03bf1a27-10bc-4362-aa76-4d24befeea4d
begin
	plot(log10.(trainingcurveDensity),
		minorgrid=true, title="Assembly interconnection density training curve",
		ylabel="interconnection density log10",
		xlabel="t", label=:none
	)
	vline!([30,30],linecolor=:black, opacity=0.3, linestyle=:dash, label=:none)
end

# ╔═╡ 4426b294-1313-46df-9807-eb22e903151f
md"By this measure, the assembly keeps becoming denser with more iterations."

# ╔═╡ 650281c0-bae5-441e-be11-72a22f00e19a
md"""
## Learning by Surprise

A central element of HTM learning is _surprise_.
In the following notebooks we'll present more elements of HTM's dynamics, but we can build on the previous mention of "minicolumns" to define a measure of how surprised the Region is with the input it receives at each point in time.
This measure is the number of local circuits that are `bursting`, ie. have all their neurons excited instead of just a few.

The distribution of excited neurons across local circuits in HTM carries a lot of information on how expected or ambiguous the input is.
Surprising input leads to many minicolumns bursting.
Ambiguous input leads to multiple neurons in many minicolumns firing, instead of just 1.
The fraction of bursting minicolumns to active minicolumns can be a measure of surprise.
"""

# ╔═╡ e759fdf6-63db-43f6-931b-5d38085d8fe8
bursting(R,c)= HierarchicalTemporalMemory.tm_activate(R.tm, c, R.tm.previous.Π)[2]

# ╔═╡ 9b5adcab-9d9e-4b4f-af9d-d44ed2248f29
surprise(R,x)= @chain x begin
	R.sp
    count(bursting(R,_)) / count(_)
end

# ╔═╡ f03586a0-7c45-4122-bb72-79bac82c1f13
md"""
We can monitor how surprised the region is as it follows a stimulus schedule that varies between 2 main stimuli:
"""

# ╔═╡ 63b3f98f-d8a6-4bb1-b3f7-0ea06e6769c4
x2= @chain bitrand(Nin) A(_).active;

# ╔═╡ 786f8a60-f4e7-42f5-aff1-f5a57049b2c7
flatCollect(x)= x|> Iterators.flatten|> collect

# ╔═╡ d73ae1ad-b35f-4c6b-a668-05700e12ec2b
stimulus= [
	[x for t=1:T],
	[x2 for t=1:T],
	[x.|x2 for t=1:T],
	[
		[[x,x2] for t=1:T/4]|> flatCollect,
		[x2,x2,x2,x2,x2],
		[[x,x2] for t=T/4+4:T]|> flatCollect,
	]|> flatCollect,
	[x for t=1:10],
	[x2 for t=1:10],
]|> flatCollect;

# ╔═╡ a43e8939-4e57-4307-a0f5-6a89b276df43
begin
	train_surprise!(R,x,x⁻)= begin
	    step!(R,x).active
		surprise(R,x)
	end
	Random.seed!(0)
	reset!(B)
	trainingcurveSurprise= [train_surprise!(B,stimulus[t],stimulus[max(t-1,1)]) for t= 1:length(stimulus)]
	md"This cell produces the training curve."
end

# ╔═╡ 32e26605-df43-4d4a-a66e-1acc33b5da6d
begin
	plot(first.(trainingcurveSurprise),
		minorgrid=true, title="Surprise as input varies",
		ylabel="surprise",
		xlabel="t", label=:none, rightmargin=12mm
	)
	textp(x,t)= (x,.98, text(t,9,:grey))
	annotate!([
		textp(.5T, "x*"),
		textp(1.5T,"x₂*"),
		textp(2.5T,"(x.|x₂)*"),
		textp(4T, "(x,x₂)*"),
		#textp(3.5T+30, "5x,5x₂"),
	])
	sep!(x)= vline!([x,x], linecolor=:grey, linestyle=:dash, linewidth=.5, label=:none)
	sep!(T); sep!(2T); sep!(3T); sep!(5T)
end

# ╔═╡ 072d080f-8e3d-4a88-a45d-ad3731cdf97a
md"""
A full discussion of what the surprise curves show, especially for more experiments with sequences, are beyond the scope of the basic operations of assembly calculus.
The important note though is that the region settles into a comfortable anticipation of both inputs `x`, `x₂`.
"""

# ╔═╡ Cell order:
# ╠═0d3bf5f6-1171-11ec-0fee-c73bb459dc3d
# ╟─1bb0fcfc-2d7a-4634-9c93-263050c56a55
# ╠═2e40a09d-3220-4424-851d-8470dfa235d1
# ╟─e109962a-dfaa-438d-aa86-3e074e57f51d
# ╠═70471bdc-660e-442e-b92a-f486abd120a7
# ╟─7620c202-c6cf-44db-9cda-13642e28b45a
# ╟─043f8da8-f17a-4486-9c49-7d5ecc75457f
# ╠═bbaf0a64-2471-4106-a27c-9dd5a4f58c7c
# ╠═9c7c7c2b-b020-4a83-8917-909470331be7
# ╟─334b271d-247c-465c-ae82-f91007a6d9d0
# ╠═4833f12f-3eac-407c-a9ce-0d1aea216077
# ╠═45d8e50b-8baa-4da2-8cd0-5c81d634c450
# ╟─91fa8319-1702-4cad-8ae8-0c1ed51a7bb8
# ╠═3ee589ed-0fe0-4328-88e8-d2fb6a98bf28
# ╟─6790a86a-3990-4c6a-878f-c18779f8d48d
# ╠═5458c6ae-ace2-4928-8a9b-ef919a1e97cb
# ╟─6b56422c-5f9b-41d8-a080-ef5ca3d7db7b
# ╟─3229da2a-92f3-4064-bb80-cbd9f5523d7d
# ╟─71f21302-04db-41aa-a259-f3be2b7bc271
# ╟─5e0ac810-1689-4bfb-88a8-85504cd821b6
# ╠═1015f629-9818-4dbb-b574-97c552e96164
# ╟─62d41be1-2970-48e1-b689-c2ecf1ab10ce
# ╠═d08cce9d-dd9c-4530-82ae-bf1db8be3281
# ╟─6a20715a-9dc9-461a-b6a2-f02522e277f0
# ╠═18e6d1c5-cd7a-4c44-a186-3ec2404f4ffc
# ╠═d59c79f9-c570-4f71-a379-6808bc189e02
# ╠═e63a5860-587c-4aff-91be-a894b3443943
# ╠═c866c71e-5dce-4af1-84bb-045815d1acb9
# ╠═75560584-1631-4c93-8808-e269d345e1c8
# ╠═456021aa-3dbf-4102-932f-43109905f9eb
# ╟─56e0da0b-8858-459a-9aae-0eba4797cc06
# ╠═cf83549c-03b1-4f7d-be29-d3a29da3e8f7
# ╠═85dfe3ed-cf70-4fde-a94a-f70c3fe4abf6
# ╠═71a2a9e5-f5cd-4c61-bfc7-4ee950897e16
# ╟─2b6d959f-63bb-4d42-a34d-70d93f6a1636
# ╠═3a5e3aa0-e34e-457c-a405-9bc95cd88910
# ╠═5ec611f4-d3ea-4c5d-9221-4dcf2da3e18c
# ╟─3df40867-674b-405a-a05c-0733a4caaa67
# ╟─c3896446-3a24-413a-afab-6a67b07092f5
# ╠═dcbe4617-7314-42c1-87bf-57ca82554a20
# ╟─059f3e13-56f6-4dc6-9570-19c74996d1ef
# ╠═8e9003fa-c822-4508-9862-58a816d9242d
# ╠═d4cfbfe7-ba6f-4936-8854-8d49277657d2
# ╟─dfff3823-d750-40e1-9ddf-769cb2c7e575
# ╠═ed61ded0-d665-46a9-bed8-cd83cb0a545a
# ╟─5bd5049c-de0d-4838-9c82-9cef604e650e
# ╟─b886095d-d449-4334-99cf-44b800bb4fc4
# ╟─8fdacada-55df-4abd-9501-7405e608529b
# ╟─6a5298f1-dc8a-4422-b034-70df607c0456
# ╟─16fa7cce-a1db-4ffd-b769-3f0109305590
# ╟─e9a82f4d-1866-45a2-a7e0-a5e08ce5e0fb
# ╠═03bf1a27-10bc-4362-aa76-4d24befeea4d
# ╟─4426b294-1313-46df-9807-eb22e903151f
# ╟─650281c0-bae5-441e-be11-72a22f00e19a
# ╟─e759fdf6-63db-43f6-931b-5d38085d8fe8
# ╟─9b5adcab-9d9e-4b4f-af9d-d44ed2248f29
# ╟─f03586a0-7c45-4122-bb72-79bac82c1f13
# ╠═63b3f98f-d8a6-4bb1-b3f7-0ea06e6769c4
# ╠═786f8a60-f4e7-42f5-aff1-f5a57049b2c7
# ╠═d73ae1ad-b35f-4c6b-a668-05700e12ec2b
# ╟─a43e8939-4e57-4307-a0f5-6a89b276df43
# ╠═32e26605-df43-4d4a-a66e-1acc33b5da6d
# ╠═072d080f-8e3d-4a88-a45d-ad3731cdf97a
