# =====================================================================================
#  Why "if-else" cannot go straight into a solver
#
#  The smallest possible instance of the problem this tutorial is about. No network, no
#  OPF, no time series: one inverter, one voltage, one reactive output, two equations.
#
#      v = V0 + K q            the grid   — the terminal voltage the inverter sees moves
#                                           with its own reactive injection
#      q = q(v)                the inverter — the IEEE 1547 Volt-VAr curve
#
#  Both are two lines of arithmetic, and the pair has exactly one solution. Getting a
#  solver to find it is where the trouble starts, and this script walks through it:
#
#      1  write the curve as a Julia `if`-`else` and hand it to JuMP   -> does not build
#      2  dodge that by evaluating the `if` on the previous iterate    -> never converges
#      3  write the curve as algebra the solver can see                -> solved, exactly
#
#  Only steps 1 and 2 are about the *droop*. They are about `if`, and the same two walls
#  stand in front of Big-M, Lambda/SOS2 and Heaviside alike, which is why all three exist.
#
#  Tutorial:
#  https://ra-emami.github.io/SmartInverterDOPF.jl/dev/tutorial_voltvar/
#
#  Run:   julia --project=. ifelse_fails.jl
# =====================================================================================

using JuMP, HiGHS, Printf

# ───────────────────────────────────────────────────────────────── the inverter's curve ──
const VBP    = [0.88, 0.90, 0.97, 1.00, 1.02, 1.10]   # IEEE 1547 breakpoints, p.u.
const QSHAPE = [1.0, 1.0, 0.0, 0.0, -1.0, -1.0]       # q / q̄ at each breakpoint
const QBAR   = 1.0                                    # reactive capability, p.u.

const A1 = -QBAR / (VBP[3] - VBP[2])    # slope of the first sloped segment
const A2 = -QBAR / (VBP[5] - VBP[4])    # slope of the second, and the steep one: -50

# ─────────────────────────────────────────────────────────────────────────── the grid ──
# A single sensitivity standing in for the whole network: absorbing VArs pulls the
# terminal voltage down, injecting them pushes it up.
const V0 = 1.035    # p.u. terminal voltage at q = 0 — a sunny afternoon, voltage already high
const K  = 0.03     # p.u. voltage per p.u. of reactive injection

# The curve, written the way anyone would write it first.
function q_droop(v)
    if     v <= VBP[2]
        return  QBAR                    # saturated, full injection
    elseif v <= VBP[3]
        return  A1 * (v - VBP[3])       # sloping down to zero
    elseif v <= VBP[4]
        return  0.0                     # dead-band
    elseif v <= VBP[5]
        return  A2 * (v - VBP[4])       # sloping down again — the steep segment
    else
        return -QBAR                    # saturated, full absorption
    end
end

banner(s) = (println(); println("=" ^ 86); println("  ", s); println("=" ^ 86))

# ═══════════════════════════════════════════ 1. hand the `if`-`else` straight to JuMP ══
banner("1  The obvious thing: put q_droop(v) in the model")

model = Model(HiGHS.Optimizer)
@variable(model, 0.90 <= v <= 1.10)
@variable(model, -QBAR <= q <= QBAR)
@constraint(model, v == V0 + K * q)          # the grid: fine, it is linear

try
    @constraint(model, q == q_droop(v))      # the droop: this is the line that fails
    println("  built without error — unexpected!")
catch err
    # JuMP's message is a page long; the first line is the whole story.
    println("  @constraint(model, q == q_droop(v))  threw:\n")
    println("    ", first(split(sprint(showerror, err), '\n')))
    println("""

  `v` is a decision variable, not a number. The solver has not chosen it yet, so
  `v <= 0.90` has no truth value to branch on, and the body of the `if` never runs.
  Branching is not an algebraic constraint, and rewriting the `if` does not change that.
  (JuMP goes on to say so at length: "You cannot write a model like this. You must
  formulate your problem as a single optimization problem.")
""")
end

# ══════════════════════════════════════════════ 2. the workaround everyone tries next ══
banner("2  The workaround: evaluate the `if` at the previous voltage, and repeat")

println("  iterate    v          q")
v_it = V0
for k in 1:12
    global v_it
    q_it = q_droop(v_it)                 # the `if` runs on a number, so this is legal
    v_it = V0 + K * q_it                 # ... and the grid answers back
    @printf("     %2d    %.4f    %+.4f\n", k, v_it, q_it)
end
@printf("""

  It never settles. The voltage lands above the last breakpoint, the inverter saturates,
  the voltage is pulled below it, the inverter backs off, and the voltage goes straight
  back up. The loop is stuck in a two-cycle because the active segment is steep enough
  that |K x slope| = %.1f > 1 — and a steep droop is exactly the case worth solving.
""", abs(K * A2))

# ═══════════════════════════════════════════════ 3. the same curve, written as algebra ══
banner("3  The fix: give the solver the whole curve at once (Lambda / SOS2)")

m = Model(HiGHS.Optimizer)
set_silent(m)
@variable(m, 0.90 <= v <= 1.10)
@variable(m, -QBAR <= q <= QBAR)
@variable(m, 0 <= lam[1:6] <= 1)             # one weight per breakpoint
@variable(m, z[1:5], Bin)                    # one binary per segment

@constraint(m, sum(lam) == 1)                # the weights are a convex combination
@constraint(m, sum(z) == 1)                  # exactly one segment is active
@constraint(m, lam[1] <= z[1])               # and only its two ends may carry weight
@constraint(m, [b in 2:5], lam[b] <= z[b-1] + z[b])
@constraint(m, lam[6] <= z[5])

@constraint(m, v == sum(lam[b] * VBP[b] for b in 1:6))            # the point on the curve
@constraint(m, q == QBAR * sum(lam[b] * QSHAPE[b] for b in 1:6))
@constraint(m, v == V0 + K * q)                                   # and the grid

optimize!(m)

vs, qs = value(v), value(q)
@printf("\n  status        %s\n", termination_status(m))
@printf("  v             %.4f p.u.\n", vs)
@printf("  q             %+.4f p.u.\n", qs)
@printf("  on the curve? q_droop(v) = %+.4f   (residual %.2e)\n",
        q_droop(vs), abs(qs - q_droop(vs)))
println("""

  One solve, no iteration, and the answer sits on the IEEE 1547 curve to machine
  precision. Nothing about the curve changed between step 1 and step 3 — only whether
  the solver could see all of it at once.
""")
