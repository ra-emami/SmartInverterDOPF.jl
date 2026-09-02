"""
    SmartInverterDOPF

Modeling smart inverters in distribution optimal power flow, by embedding the
IEEE 1547 Volt-VAr (Q-V) droop curve of each inverter into the OPF itself.

The package provides two distribution OPF host models and three interchangeable
encodings of the droop law, chosen independently of one another.

Droop encodings (`method`):

| method | class | requires |
|:--|:--|:--|
| `:bigm` | mixed-integer linear | an MILP solver |
| `:lambda` | mixed-integer linear (SOS2) | an MILP solver |
| `:heaviside` | non-smooth nonlinear | an NLP solver |

All three describe the same curve and return the same dispatch; they differ in the
solver technology they need and in how they scale.

Host models (`host`):

| host | model | accuracy | solve |
|:--|:--|:--|:--|
| `:ivacopf` (default) | current-voltage AC-OPF | very accurate, near-exact AC | iterative successive linearisation |
| `:lindistflow` | linearised branch flow | approximate: losses dropped | run once, much faster |

```julia
using SmartInverterDOPF, Gurobi
case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)                      # IVACOPF
fast = solve_dopf(case, Gurobi.Optimizer; method = :lambda, host = :lindistflow) # LinDistFlow
```

!!! warning
    HiGHS does not complete the successive-linearisation loop on the bundled case: it
    solves the first pass, then reports `INFEASIBLE` on the second, for both MILP
    encodings, at every MIP gap tried. Gurobi solves the same sequence to proven
    optimality.
"""
module SmartInverterDOPF

using JuMP
using JSON3

export Case, DroopCurve, DOPFResult,
       load_case, solve_dopf, base_case_voltages,
       ieee1547_curve, droop_q, kWh, nbus, ndg

include("case.jl")
include("droop.jl")
include("dopf.jl")

end # module
