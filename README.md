# SmartInverterDOPF.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Modeling Smart Inverters in Distribution Optimal Power Flow

A smart inverter does not accept a reactive-power set-point. It follows a Volt-VAr
curve, deciding from its own terminal voltage how much reactive power to inject or
absorb. An OPF that ignores that curve returns a dispatch the inverter is never going to
deliver. This package puts the curve inside the optimisation, three different ways, and
runs each of them against four distribution OPF host models: two single-phase, two
three-phase, linear and near-exact in each pair. The encodings agree; the hosts do not,
and an exact AC power flow decides between them.

**[Read the tutorial →](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/tutorial_voltvar/)**

## Quick start

```julia
using Pkg
Pkg.add(url = "https://github.com/ra-emami/SmartInverterDOPF.jl")
Pkg.add(["Gurobi", "Ipopt"])                    # solvers
Pkg.add(["JuMP", "JSON3", "Plots", "Printf"])   # modelling, case files, figures, tables
```

Every package a `using` line names has to be installed into the active environment
first, not only the solvers. The tutorial's
[Prerequisites](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/tutorial_voltvar/#Prerequisites)
section lists the full set.

```julia
using SmartInverterDOPF, Gurobi

case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)   # host = :ivacopf (default)

kWh(case, sum(res.PVC))     # PV energy curtailed over the day, kWh
extrema(res.V)              # voltage range across the feeder, p.u.
```

## The single-phase host models

`method` picks how the droop curve is encoded; `host` picks the network model it sits
inside. They are independent: the droop constraints are identical in both hosts.

| `host` | model | accuracy | solve |
|:--|:--|:--|:--|
| `:ivacopf` *(default)* | current-voltage AC-OPF ([doi:10.1109/OJIA.2024.3367547](https://doi.org/10.1109/OJIA.2024.3367547), extended in [doi:10.1016/j.epsr.2026.113613](https://doi.org/10.1016/j.epsr.2026.113613)) | very accurate, near-exact AC | **iterative**: re-linearised until the residual of the exact power-flow identity clears a tolerance |
| `:lindistflow` | linearised branch flow ([doi:10.1109/SMARTGRID.2010.5622021](https://doi.org/10.1109/SMARTGRID.2010.5622021)) | approximate: losses dropped, small voltage deviations assumed | run **once**, much faster and far lower computational effort |

```julia
solve_dopf(case, Gurobi.Optimizer; method = :lambda, host = :lindistflow)
```

> **Use `:ivacopf` for quantitative work.** Audited against an exact AC power flow on the
> bundled case, the IVACOPF dispatch reproduces the true solution to ~1e-9 p.u. and sits
> on the droop curve to ~1e-7. The LinDistFlow dispatch is off the droop by 6% of
> inverter rating and puts 17 of the 96 time steps below the lower voltage limit, so it
> is not deliverable. `:lindistflow` is for a fast first look and for warm-starting.

The cheap host also makes the accurate one cheaper. `warm_start = :lindistflow` seeds
the IVACOPF linearisation from the LinDistFlow solution instead of a flat profile, which
halves the passes and is ~1.7× faster overall, for the same answer:

```julia
solve_dopf(case, Gurobi.Optimizer; method = :lambda, warm_start = :lindistflow)
```

> **Solver note.** The committed results are generated with **Gurobi**, which is
> [free for academic users](https://www.gurobi.com/academia/academic-program-and-licenses/).
> We also tried the open-source MILP solvers HiGHS and GLPK; neither worked out. One
> returned an infeasible status inside the successive-linearisation loop, the other was
> too slow to finish. If no MILP licence is available, the `:heaviside` encoding needs
> only Ipopt, which is open source, and reaches the same answer.

## The three encodings

The Volt-VAr law is a five-segment piecewise-linear function, a definition by cases,
which is exactly what a solver cannot read. Three standard rewrites make it tractable:

| `method` | class | idea | needs |
|:--|:--|:--|:--|
| `:bigm` | MILP | one binary per segment activates that segment's voltage window and affine law | MILP solver |
| `:lambda` | MILP | operating point as a convex combination of breakpoints, SOS2 forcing adjacency | MILP solver |
| `:heaviside` | NLP | segment masks from unit steps, summed into one closed-form expression | NLP solver |

All three are exact and, on the bundled case study, return the same dispatch to within
solver tolerance. They differ in the solver technology they demand and in how they scale
with inverters × time steps.

## Case studies

Three feeders, all over 24 h at 15-minute resolution (96 steps), all minimising PV
curtailment, all with voltages held inside `[0.95, 1.05]` p.u.

| feeder | phases | size | inverters | used for |
|:--|:--|:--|:--|:--|
| IEEE 33-bus | single | 33 buses | 3, at buses 7 / 18 / 33 | the main case study; the package's `load_case()` |
| `network_5_Feeder_2` | three | 194 buses, 18 single-phase loads split 4/5/9 | 12 in four size classes, 84 kW | the three-phase case study |
| `network_17_Feeder_6` | three | 3856 buses, 223 single-phase loads | 12 | the scalability check |

Both three-phase feeders are real Electricity North West low-voltage networks, Kron-reduced
to three wires; sources and licence are in
[`examples/three_phase/README.md`](examples/three_phase/README.md).

Headline results, all three encodings agreeing to solver tolerance within each host:

| case | host | curtailed | droop residual at the **true** AC voltage |
|:--|:--|--:|--:|
| 33-bus | IVACOPF | 3022.0 kWh | 7.9e-08 |
| 33-bus | LinDistFlow | 60.5 kWh | 6.1e-02, **not deliverable** |
| 194-bus | IVACOPF | 42.69 kWh | 2.8e-11 |
| 194-bus | LinDist3Flow | 46.32 kWh | 8.0e-03 |

The last column is the test that separates the hosts: take each dispatch, solve the exact
AC power flow for it, and ask whether the inverters would really have produced those VArs.

## Repository layout

```
src/                    the package: case data, droop encodings, DOPF hosts
data/                   IEEE 33-bus feeder, load profiles, solar profile (JSON)
scripts/                generate_results.jl: recomputes the single-phase results
examples/single_phase/  6 standalone scripts: {ivacopf,lindistflow} × {bigm,lambda,heaviside}
examples/three_phase/   6 standalone scripts: {LinDist3Flow,IVACOPF3Ph} × {BigM,Lambda,Heaviside}
                        plus generate_results.jl, scalability.jl, plot_network.jl
docs/                   Documenter site; builds without any solver
test/                   test suite
```

Every standalone script is self-contained and shares its skeleton verbatim with its
siblings, so a `diff` between any two shows only the droop block, or only the network
model.

## Three-phase examples

[`examples/three_phase/`](examples/three_phase) runs the same three encodings on a real
unbalanced LV feeder, against **two** three-phase hosts:

| | Big-M | Lambda / SOS2 | Heaviside |
|:--|:--|:--|:--|
| **LinDist3Flow**: multiphase linearised branch flow, one pass | `LinDist3Flow_BigM.jl` | `LinDist3Flow_Lambda.jl` | `LinDist3Flow_Heaviside.jl` |
| **IVACOPF**: three-phase current-voltage AC-OPF, successive linearisation | `IVACOPF3Ph_BigM.jl` | `IVACOPF3Ph_Lambda.jl` | `IVACOPF3Ph_Heaviside.jl` |

```bash
julia --project=examples/three_phase examples/three_phase/IVACOPF3Ph_Lambda.jl
```

Feeder, horizon and fleet come from the environment (`TP_CASE`, `TP_STEPS`, `TP_NPV`, and
for IVACOPF `TP_WARMSTART`, `TP_TOL`, `TP_MAXITER`, `TP_IMAXSEG`), so the same model runs
on a different network without editing anything.

### Does it scale?

`scalability.jl` reruns the LinDist3Flow scripts on the 3856-bus feeder. The mixed-integer
encodings carry a 3.3-million-variable model over the full day in about a minute, and the
binary count does not move between feeders: it depends on inverters × time steps, not on
network size. The integer-free encoding is the one that breaks: Ipopt gives up at the full
horizon and needs a shortened day to finish. All three stay exact wherever they finish.

## Regenerating the documentation results

The documentation is built from precomputed results committed under
`docs/src/assets/results/`, so building the docs needs no optimisation solver. To
recompute them:

```bash
julia --project=scripts scripts/generate_results.jl
```

and, for the three-phase section, which runs all six scripts, both hosts:

```bash
julia --project=examples/three_phase examples/three_phase/generate_results.jl
```

`TP_HOSTS=ivacopf` or `TP_HOSTS=lindist3flow` regenerates one family only.

## Building the documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

MIT. See [LICENSE](LICENSE).
