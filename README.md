# Marginal Treatment Effects

Replication package for:

1. "Using Instrumental Variables for Inference About Policy Relevant Treatment Parameters"
    Mogstad, Santos, and Torgovitsky (2018, Econometrica)
2. "Identification and Extrapolation of Causal Effects with Instrumental Variables"
    Mogstad and Torgovitsky (2018, Annual Review of Economics)

To load/install dependencies:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate() # download the correct packages
```

Then run with

```julia
using MarginalTreatmentEffects
menu("/path/to/where/you/want/to/save")
```

If you want to try to automatically call `latexpdf` to compile the TikZ figures, do this:
```julia
menu("/path/to/where/you/want/to/save", compile = true)
```

## Disclaimer

The solution to the linear programs in MST (2018) and MT (2018) need not be
unique. While the bounds on the target parameters shouldn't change, the MTRs
that achieve these bounds may be different when running the code on different
machines.
