# Why "if-else" cannot go straight into a solver

The smallest possible instance of the problem the [tutorial](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/tutorial_voltvar/)
is about. No network, no OPF, no time series: one inverter, one voltage, one reactive
output, two equations.

```math
v = V_0 + K q          # the grid: the terminal voltage moves with the inverter's own VArs
q = q(v)               # the inverter: the IEEE 1547 Volt-VAr curve
```

Both are two lines of arithmetic and the pair has exactly one solution,
`v = 1.0140`, `q = -0.7000`. Getting a solver to find it is where the trouble starts.

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. ifelse_fails.jl
```

Needs only JuMP and HiGHS, both open source. No licence, nothing to download beyond the
two packages.

## What the script shows

**1. Write the curve as a Julia `if`-`else` and hand it to JuMP.** The model does not
even build:

```
Cannot evaluate `<=` between a variable and a number.
```

`v` is a decision variable. The solver has not chosen it yet, so `v <= 0.90` has no truth
value to branch on and the body of the `if` never runs. Branching is not an algebraic
constraint, and rewriting the `if` does not change that.

**2. Dodge it by evaluating the `if` at the previous voltage, and repeating.** Now the
`if` runs on a number, so it is legal — and it never converges:

```
  iterate    v          q
      1    1.0050    -1.0000
      2    1.0275    -0.2500
      3    1.0050    -1.0000
      4    1.0275    -0.2500
      ...
```

The voltage lands above the last breakpoint, the inverter saturates, the voltage is
pulled below it, the inverter backs off, and the voltage goes straight back up. The loop
is stuck in a two-cycle because the active segment is steep enough that
`|K x slope| = 1.5 > 1` — and a steep droop is exactly the case worth solving.

**3. Write the same curve as algebra the solver can see.** The Lambda/SOS2 encoding, in
nine lines of JuMP:

```
  status        OPTIMAL
  v             1.0140 p.u.
  q             -0.7000 p.u.
  on the curve? q_droop(v) = -0.7000   (residual 2.22e-15)
```

One solve, no iteration, and the answer sits on the curve to machine precision. Nothing
about the curve changed between step 1 and step 3 — only whether the solver could see all
of it at once.

## Where this goes next

Steps 1 and 2 are not really about the droop; they are about `if`. The same two walls
stand in front of every encoding, which is why the tutorial carries three of them:

| encoding | idea | class |
|:--|:--|:--|
| **Big-M** | a binary per segment activates that segment's law and voltage range | MILP |
| **Lambda / SOS2** | the operating point is a blend of two adjacent breakpoints — used above | MILP |
| **Heaviside** | segment masks built from unit steps, no integers at all | NLP |

All three are exact, all three return the same dispatch, and all three drop unchanged
into a real distribution OPF: see [`examples/single_phase`](../single_phase) and
[`examples/three_phase`](../three_phase).
