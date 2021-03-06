### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 201a9fd8-7f37-48ea-b8ad-0cb1ac196c43
using HierarchicalTemporalMemory

# ╔═╡ bcee52a7-7156-4b8e-8000-f6ef2d7a7078
using Setfield

# ╔═╡ 8429596b-8765-4548-8300-b5449889395d
using Chain

# ╔═╡ 2d930a14-b39e-462b-8e4d-b297573e1ea3
using Statistics

# ╔═╡ 0ada6b68-73d3-11ec-33d7-3b010e6c893e
md"""
# Assembly operations library

This notebook is a library to be imported by the other notebooks, accumulating pieces of code developed there so that they don't have to be re-written.
It focuses on the assembly operations.
"""

# ╔═╡ ea459f12-2b04-41a8-946f-a3ac3d9d040c
"""
    reldistance(yᵢ,yⱼ)

Relative distance between 2 sparse binary vectors bounded in ``[0,1]``:
``yᵢ = yⱼ ⇔ reldistance(yᵢ,yⱼ) = 0``
and
``yᵢ ⊥ yⱼ ⇔ reldistance(yᵢ,yⱼ) = 1``
"""
reldistance(yᵢ,yⱼ)= 1 - count(yᵢ .& yⱼ) * mean(inv, [count(yᵢ), count(yⱼ)])

# ╔═╡ fdbd6b4d-dad6-4b2c-b9cc-208ac34255b1
"`Sₜ(yₜ,τ)` is the moving average of neurons that have been activated in the region over the last τ timesteps, given the activations of the region over time yₜ."
Sₜ(yₜ,τ)= [ reduce((a,b) -> a.|b , @view yₜ[max(1,t-τ):t]) for t=1:length(yₜ) ]

# ╔═╡ f97ab6f6-cf5f-4b50-ac67-33da6abc3c41
"`Δrel(x)` is the relative distance between 2 successive elements of x."
Δrel(x)= @views [0; map(reldistance, x[1:end-1], x[2:end])]

# ╔═╡ 36b03905-56a6-4c31-b09d-d5d4c05be704
convergence(yₜ; τ=5)= begin
	s= Sₜ(yₜ,τ)
	Δrel(s)
end

# ╔═╡ 4c15ea56-c8ec-4d3f-bdb6-ab36cf357b7e
"""
    project!(R,x; time_to_convergence=30)

Create a projection of feedforward activation `x` on region `R`.
The region `R` adapts to the stimulation.

See also [01-intro_projection.jl](./open?path=notebooks/01-intro_projection.jl)
"""
project!(R,x; time_to_convergence=30)= Iterators.drop((step!(R,x).active for t= 1:time_to_convergence), time_to_convergence-1)|> first

# ╔═╡ 37056ba5-93a4-49f7-b285-9bf750ea52ae
"""
`probe(R,x)` shows the input `x` 2 times to the region.
This temporarily changes the region's predictive state, therefore more faithfully representing the region's response to a stimulus in the context of assemblies, compared to `R(x)`.
It doesn't modify `R`.
"""
probe(R,x)= begin
	_R= deepcopy(R)
	step!(_R,x, learn= false)
	_R(x)
end

# ╔═╡ 1598b575-7bd2-48b1-97ab-e1a41a8d1680
"`interconnectionMeasure(x,R)` counts the number of distal synapses between neurons of the given activation pattern `x` in region `R`."
interconnectionMeasure(x, R)= x' * distalSynapses(R) * x

# ╔═╡ 779253f0-0526-4ea4-8eb3-5aa0de55c14c
"`interconnectionDensity(x,R)` is the ratio of connected synapses whose pre- and post-synaptic neurons are inside the activation set versus those whose pre- synaptic neurons are outside it."
interconnectionDensity(x,R)= x' * distalSynapses(R) * x / ((.!x)' * distalSynapses(R) * x)

# ╔═╡ 08fd02b1-c020-4b7f-86c7-bcfd5c9d1a40
"`subset(x,p)` activates a random subset of neurons with ratio `p` from activation `x`. Note: `0≤p≤1`"
subset(x,p)= @chain [i for i= HierarchicalTemporalMemory.Truesof(x)] begin
	shuffle(_)[1 : Int(p*count(x))]
	HierarchicalTemporalMemory.bitarray(_,length(a))
end

# ╔═╡ 142ecb6b-01d0-4c30-9ad8-60f7f697e173
"`sparseBitrand(dims, sparsity)` like `bitrand` creates a random bitvector, but with only a small percentage `sparsity` of `Trues`."
sparseBitrand(dims, sparsity)= HierarchicalTemporalMemory.bitarray(randsubseq(1:dims, sparsity), dims)

# ╔═╡ c199f632-b04e-4a29-acf1-de5ae63429fd
"`overlap(x,y)` of 2 assemblies in the same region is the number of neurons that are active in both. Symbol: `x&y`"
overlap(x,y)= count(x .& y)

# ╔═╡ 6bd2e07f-3ea7-4e30-b666-7f0c153e5ee3
"`minicolumnOverlap(a,b,R)` is the number of minicolumns that are active in both a,b (of region R)"
minicolumnOverlap(a,b,R)= overlap(
	any(reshape(a,R.tm.params.k,:),dims=1),
	any(reshape(b,R.tm.params.k,:),dims=1)
)

# ╔═╡ 29d5f907-a110-490c-94a0-fe707a96b99d
"""
    bursting(R,c)

For each minicolumn, show if it is bursting (surprise) in response to minicolumn activation `c`.

### Example

Create an assembly and get the region's bursting minicolumns when presented with the same input again.
```
r= Region(params_input.sp, params_input.tm)
aᵢₙ= bitrand(Nin)
assembly= project!(r, aᵢₙ)
bursting(r, r.sp(aᵢₙ))
```
"""
bursting(R,c)= HierarchicalTemporalMemory.tm_activate(R.tm, c, R.tm.previous.Π)[2]

# ╔═╡ 8f57d50e-dccc-4440-b252-a7f5ef34d084
"`surprise(R,x)` counts the fraction of [`bursting`](@ref) minicolumns in response to region input `x`"
surprise(R,x)= @chain x begin
	R.sp
    count(bursting(R,_)) / count(_)
end

# ╔═╡ ee8b169e-a096-40fb-99c8-0365cdbed953
flatCollect(x)= x|> Iterators.flatten|> collect

# ╔═╡ 11d7cda5-4103-4a94-b7a7-1fbfbaaa8c5c
md"""
## HTM regions

This section defines standard HTM regions that were tuned to work well with assembly calculus operations.

First, we define parameters that have been tuned to work well with assembly calculus.
"""

# ╔═╡ 1cde9527-e7fc-4fc3-9759-7616b40ac2ad
begin
	Nin= 1e3|> Int        		# input size
	Nn= 20e4             	    # number of neurons in each area
	k= 15                 		# neurons per minicolumn
	# calibrate how many dendritic inputs needed to fire
	# these are related to the size of patterns, which ranges between
	# `_Nc()*sparsity() * [1..k]` (from perfectly unambiguous to bursting)
	thresholds= (
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

# ╔═╡ 2698fbdb-5bd1-4902-b6d7-be1cdd9be5ad
md"""
Note: too high value for `tm_dendriteSynapses` can cause *instability*
"""

# ╔═╡ 25f3e7a1-21a6-4458-90af-349b6c3cbff4
md"Then, we produce parameter arrays mostly with these values that are used in constructing HTM Regions."

# ╔═╡ de50d108-6522-42d4-ae86-6d27b805c5d5
begin
	_Nc()= floor(Int,Nn/k) 		# number of minicolumns
	sparsity()= 1/sqrt(_Nc());
	params_input= (
		sp= SPParams(szᵢₙ=Nin, szₛₚ=_Nc(), s= sparsity(), prob_synapse=5e-3,
			p⁺_01= learnrate.prox_p⁺, p⁻_01= learnrate.prox_p⁻,
			enable_local_inhibit=false),
		tm= TMParams(Nc=_Nc(), k=k,
			p⁺_01= learnrate.dist_p⁺, p⁻_01= learnrate.dist_p⁻, LTD_p⁻_01= learnrate.dist_LTD_p⁻,
			θ_stimulus_learn=thresholds.tm_learn,
			θ_stimulus_activate=thresholds.tm_activate,
			synapseSampleSize=thresholds.tm_dendriteSynapses)
	)
	# parameters for "merge" or "work" region, where assemblies will be formed
	params_M= @set params_input.sp.szᵢₙ= params_input.tm.Nₙ
end

# ╔═╡ 6cf24775-ffcd-4c05-8219-d4850d21cefa
_Nc() .* (1,sparsity())

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
HierarchicalTemporalMemory = "db654dba-670c-566a-9226-64313b881ded"
Setfield = "efcf1570-3423-57d1-acb7-fd33fddbac46"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Chain = "~0.4.10"
HierarchicalTemporalMemory = "~0.3.0"
Setfield = "~0.8.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "926870acb6cbcf029396f2f2de030282b6bc1941"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.4"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.HierarchicalTemporalMemory]]
deps = ["Chain", "ImageFiltering", "IterTools", "Lazy", "LinearAlgebra", "Parameters", "Random", "Setfield", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "a2b935a3a4084d4f4e0964587879089519fcea94"
uuid = "db654dba-670c-566a-9226-64313b881ded"
version = "0.3.0"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2ae4fe21e97cd13efd857462c1869b73c9f61be3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═201a9fd8-7f37-48ea-b8ad-0cb1ac196c43
# ╠═bcee52a7-7156-4b8e-8000-f6ef2d7a7078
# ╠═8429596b-8765-4548-8300-b5449889395d
# ╠═2d930a14-b39e-462b-8e4d-b297573e1ea3
# ╟─0ada6b68-73d3-11ec-33d7-3b010e6c893e
# ╠═ea459f12-2b04-41a8-946f-a3ac3d9d040c
# ╠═fdbd6b4d-dad6-4b2c-b9cc-208ac34255b1
# ╠═f97ab6f6-cf5f-4b50-ac67-33da6abc3c41
# ╠═36b03905-56a6-4c31-b09d-d5d4c05be704
# ╠═4c15ea56-c8ec-4d3f-bdb6-ab36cf357b7e
# ╠═37056ba5-93a4-49f7-b285-9bf750ea52ae
# ╠═1598b575-7bd2-48b1-97ab-e1a41a8d1680
# ╠═779253f0-0526-4ea4-8eb3-5aa0de55c14c
# ╠═08fd02b1-c020-4b7f-86c7-bcfd5c9d1a40
# ╠═142ecb6b-01d0-4c30-9ad8-60f7f697e173
# ╠═c199f632-b04e-4a29-acf1-de5ae63429fd
# ╠═6bd2e07f-3ea7-4e30-b666-7f0c153e5ee3
# ╟─29d5f907-a110-490c-94a0-fe707a96b99d
# ╠═8f57d50e-dccc-4440-b252-a7f5ef34d084
# ╠═ee8b169e-a096-40fb-99c8-0365cdbed953
# ╟─11d7cda5-4103-4a94-b7a7-1fbfbaaa8c5c
# ╠═1cde9527-e7fc-4fc3-9759-7616b40ac2ad
# ╟─2698fbdb-5bd1-4902-b6d7-be1cdd9be5ad
# ╟─25f3e7a1-21a6-4458-90af-349b6c3cbff4
# ╠═de50d108-6522-42d4-ae86-6d27b805c5d5
# ╠═6cf24775-ffcd-4c05-8219-d4850d21cefa
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
