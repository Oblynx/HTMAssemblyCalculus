using DrWatson; @quickactivate
using Revise
import HierarchicalTemporalMemory as HTM
using Random, Chain, Plots

# Region 1: feedforward + recurrent
Nin= 1e3|> Int
Nn= 4e5
k= 10
Nc()= floor(Int,Nn/k)
sparsity()= 1/sqrt(Nc())
Nₙ()= tm.params.Nₙ
sp= HTM.SpatialPooler(HTM.SPParams(szᵢₙ=Nin, szₛₚ=Nc(), s= sparsity(), prob_synapse=1e-3, enable_local_inhibit=false))
tm= HTM.TemporalMemory(HTM.TMParams(Nc=Nc(), k=k, p⁺_01= .11, θ_stimulus_learn=15, θ_stimulus_activate= 22))

# Region 2
sp2= HTM.SpatialPooler(HTM.SPParams(szᵢₙ=Nₙ(), szₛₚ=Nc(), s= sparsity(), prob_synapse=1e-3, enable_local_inhibit=false))
tm2= HTM.TemporalMemory(HTM.TMParams(Nc=Nc(), k=k, p⁺_01= .11, θ_stimulus_learn=15, θ_stimulus_activate= 22))

## Send SDRs from 1 -> 2
T= 35; experiments=5
a= Matrix(undef,T,experiments)
for e= 1:experiments
  input= bitrand(Nin)
  z= @chain input sp tm (_).active
  # Activation of first region as input of the 2nd
  for t= 1:T
    t%10==0 && @info "t=$t"
    # Column activation + feedforward plasticity
    c= HTM.step!(sp2,z)
    a[t,e]= HTM.step!(tm2,c).active
  end
end

similarity= [ 2count(a[i,e] .& a[i+1,e])/(count(a[i,e]) + count(a[i+1,e]))
    for i=4:T-1, e=1:experiments
]
plot(similarity, minorgrid=true)
