# Marginal Treatment Effects

Replication package for:

1. "Using Instrumental Variables for Inference About Policy Relevant Treatment Parameters"
    Mogstad, Santos, and Torgovitsky (2018, _Econometrica_)
2. "Identification and Extrapolation of Causal Effects with Instrumental Variables"
    Mogstad and Torgovitsky (2018, _Annual Review of Economics_)
3. "Policy Evaluation with Multiple Instrumental Variables"
    Mostad, Torgovitsky, and Walters (Forthcoming, _Journal of Econometrics_)

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

The solutions to the linear programs are typically not unique and can depend intricately on both software versions and hardware.
As a consequence, figures that graph optimal MTR functions (such as Figures 2--7 of MST) may not reproduce, and might reproduce differently on different machines.
The bounds on the target parameters are unique, however, and should reproduce on different machines.

## Errata

### MST (2018, Econometrica)

- Figure 4: The IV Slope entries in the legend correctly uses the support of the instruments:
  - 𝟙(Z = 2) is now 𝟙(Z = 1)
  - 𝟙(Z = 3) is now 𝟙(Z = 2)
- Figures 5, 6, 7: The entries in the legend correctly use the support of the instruments:
  - 𝟙(Z = 1) is now 𝟙(Z = 0)
  - 𝟙(Z = 2) is now 𝟙(Z = 1)
  - 𝟙(Z = 3) is now 𝟙(Z = 2)
- Figure 7: The title was changed from "Order 9 polynomial bounds" to "9th degree polynomial bounds".

### MT (2018, Annual Review of Economics)

- Figure 3: The label on the top tick of Figure 3 was changed from $p(1) - p(2)$ to $p(2) - p(1)$.
