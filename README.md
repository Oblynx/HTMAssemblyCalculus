# A calculus for brain computation based on the Hierarchical Temporal Memory

Compare the [Hierarchical Temporal Memory](https://github.com/Oblynx/HierarchicalTemporalMemory.jl) with another model of brain computation based on [neuron assemblies](https://www.pnas.org/content/117/25/14464).

**Neuron Assembly Calculus** is a computational framework for cognitive function.
It describes a dynamical system of brain areas with random connections between excitatory neurons, with Hebbian plasticity, where only the top-k most activated neurons fire.
From these properties _neuron assemblies_ emerge: clusters of highly interconnected neurons in the same area,
which can be created through programmatic operations.
The authors probabilistically prove the convergence of these operations based on the model's structure and connectome (ref: [Neuron assembly calculus], [A biologically plausible parser]).

**Hierarchical Temporal Memory** (HTM) is a biologically constrained model of brain computation.
It too is a dynamical system of brain areas with excitatory neurons, Hebbian plasticity and top-k activation.
However, the dynamics are more strongly constrained. For example,
neurons are arranged in local circuits with local competition, and have multiple dendrites, where they receive less-stimulating input ([Why neurons have thousands of synapses]).
This model of a brain area's dynamics underpins the [Thousand Brains theory] that attempts to explain human neocortical computation.

## Questions

The 2 models are based on similar, but not identical assumptions.
Contrasting these 2, a few questions emerge that can draw insight on both models:

1. Does the Hierarchical Temporal Memory create neuron assemblies? Can we manipulate them with the same operations?
1. How sensitive is neuron Assembly Calculus to its precise set of assumptions? Will it hold up to the stronger biological constraints imposed by HTM?

A positive answer to these questions can strengthen both models.

# Implement neuron assemblies on HTM

This project uses the Julia implementation of HTM ([HTM.jl]) to experiment with neuron assemblies.
We implement the basic operations of Assembly Calculus to explore the conditions of their convergence by simulating the HTM dynamics.

## Notebooks

We experiment in [Pluto notebooks], which provide a reproducible environment.
The rendered HTML version is available at https://oblynx.github.io/HTMAssemblyCalculus

The notebooks so far are the following:

1. [projection.jl](notebooks/projection.jl): show the most basic operation, the projection.



[Neuron assembly calculus]: https://www.pnas.org/content/117/25/14464
[A biologically plausible parser]: https://direct.mit.edu/tacl/article/doi/10.1162/tacl_a_00432/108608/A-Biologically-Plausible-Parser
[Thousand Brains theory]: https://link.springer.com/article/10.1007/s42452-021-04715-0
[HTM.jl]: https://github.com/Oblynx/HierarchicalTemporalMemory.jl
[Why neurons have thousands of synapses]: https://www.frontiersin.org/articles/10.3389/fncir.2016.00023/full
[Pluto notebooks]: https://github.com/fonsp/Pluto.jl