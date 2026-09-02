# Current-voltage AC optimal power flow (IVACOPF), solved by successive linearisation,
# with the Volt-VAr droop plugged in as a self-contained module.

"""
    DOPFResult

Outcome of [`solve_dopf`](@ref).

| field | meaning |
|:--|:--|
| `V` | bus voltage magnitude, `[bus, hour, quarter]`, p.u. |
| `Vdg`, `Qdg`, `Pdg`, `PVC` | per-inverter voltage, reactive output, active output and curtailment, `[dg, hour, quarter]`, p.u. |
| `Ploss` | total network active loss over the day, p.u. |
| `iterations` | per-iteration solve time, objective and linearisation residual |
| `nvar`, `nbin`, `ncon` | size of the model actually handed to the solver |
| `converged` | whether the residual fell below `tol` |
| `solve_seconds` | summed solver time over all iterations |

Two fields read differently under `host = :lindistflow`, which has no outer loop:
`iterations` holds a single entry with `residual = 0.0`, and `converged` is always `true`
because there is nothing to converge: the model is linear and solved once. `Ploss` is
then a post-hoc estimate `Σ R(P² + Q²)/v²` from the solved flows rather than a modelled
quantity, since LinDistFlow drops losses from the power balance.
"""
struct DOPFResult
    V::Array{Float64,3}
    Vdg::Array{Float64,3}
    Qdg::Array{Float64,3}
    Pdg::Array{Float64,3}
    PVC::Array{Float64,3}
    Ploss::Float64
    iterations::Vector{NamedTuple{(:iter, :seconds, :objective, :residual, :status),
                                  Tuple{Int,Float64,Float64,Float64,String}}}
    nvar::Int
    nbin::Int
    ncon::Int
    converged::Bool
    solve_seconds::Float64
end

"""
    solve_dopf(case, optimizer; method = :lambda, host = :ivacopf,
               curve = ieee1547_curve(), max_iter = 15, tol = 1e-6, silent = true)

Minimise total PV curtailment over the day subject to a distribution network model and
the Volt-VAr droop encoded by `method` (`:bigm`, `:lambda` or `:heaviside`).

`host` selects the network model, and the droop encoding is identical in both:

- `:ivacopf` (default): current-voltage AC-OPF. The AC power flow enters through two
  bilinear identities, the `v·I` power balance and the `|I|²` loss, which are linearised
  about the previous iterate and refreshed until the residual of the *exact* loss
  identity falls below `tol`. Near-exact, at the cost of an outer loop.
- `:lindistflow`: the linearised branch-flow model, with losses dropped from the balance and
  the voltage drop taken as `-(R·P + X·Q)/Vnom`. Wholly linear, solved **once**, so
  `max_iter` and `tol` are ignored. Much cheaper, at the cost of accuracy. `Ploss` is
  then a post-hoc estimate (see [`DOPFResult`](@ref)).

`warm_start = :lindistflow` solves the linear host first and expands its dispatch into a
consistent complex state with an exact power-flow sweep, using that as the initial
linearisation point instead of a flat profile. On the bundled case this halves the number
of passes and is about 1.7× faster overall, and reaches the same answer.

`optimizer` must match the droop encoding: an MILP solver for `:bigm` and `:lambda`, an
NLP solver for `:heaviside`.

!!! note "Which host for quantitative work"
    `:ivacopf`. Audited against an exact AC power flow on the bundled case, its dispatch
    reproduces the true solution to ~1e-9 p.u. and sits on the droop to ~1e-7, while the
    LinDistFlow dispatch is off the droop by 6% of inverter rating and puts 17 of the 96
    time steps below the lower voltage limit. `:lindistflow` is for a fast first look and
    for warm-starting, not for reporting.

```julia
using SmartInverterDOPF, Gurobi
case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)                      # IVACOPF
warm = solve_dopf(case, Gurobi.Optimizer; method = :lambda, warm_start = :lindistflow)
fast = solve_dopf(case, Gurobi.Optimizer; method = :lambda, host = :lindistflow) # one LP/MILP
```
"""
function solve_dopf(c::Case, optimizer;
                    method::Symbol = :lambda,
                    host::Symbol = :ivacopf,
                    warm_start::Union{Nothing,Symbol} = nothing,
                    curve::DroopCurve = ieee1547_curve(),
                    max_iter::Int = 15,
                    tol::Float64 = 1e-6,
                    silent::Bool = true,
                    time_limit_sec = 3600,
                    attributes = Dict{String,Any}())

    host in (:ivacopf, :lindistflow) ||
        throw(ArgumentError("host must be :ivacopf or :lindistflow, got :$host"))
    host === :lindistflow && return _solve_lindistflow(c, optimizer; method, curve,
                                                       silent, time_limit_sec, attributes)

    HOUR_SET, QUARTER_SET = c.HOUR_SET, c.QUARTER_SET
    BUS_SET, BRANCH_SET   = c.BUS_SET, c.BRANCH_SET
    Bi_BRANCH_SET         = c.Bi_BRANCH_SET
    DG_SET, NON_DG_SET    = c.DG_SET, c.NON_DG_SET
    SLACK_SET             = [c.slack]
    R, X                  = c.R, c.X
    nb                    = length(BUS_SET)
    qbar                  = Dict(d => c.Sdg_max[d] for d in DG_SET)   # reactive capability

    # Linearisation point. A flat start (v = 1∠0, all currents zero) is far from any
    # solution, and the trajectory the successive linearisation then follows depends on
    # it. `warm_start = :lindistflow` instead solves the linear host first and expands
    # its dispatch into a consistent complex state with an exact power-flow sweep, so
    # the first linearisation is taken about a point that is already close.
    v_r_pr    = ones(nb, 24, 4)
    v_im_pr   = zeros(nb, 24, 4)
    Ibs_r_pr  = zeros(nb, 24, 4)
    Ibs_im_pr = zeros(nb, 24, 4)
    Ibr_r_pr  = Dict((br, h, q) => 0.0 for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET)
    Ibr_im_pr = Dict((br, h, q) => 0.0 for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET)
    warm_seconds = 0.0

    if warm_start === :lindistflow
        warm_seconds = @elapsed begin
            seed = _solve_lindistflow(c, optimizer; method, curve, silent,
                                      time_limit_sec, attributes)
            v_r_pr, v_im_pr, Ibs_r_pr, Ibs_im_pr, Ibr_r_pr, Ibr_im_pr =
                _sweep_state(c, seed.Pdg, seed.Qdg)
        end
    elseif warm_start !== nothing
        throw(ArgumentError("warm_start must be nothing or :lindistflow, got :$warm_start"))
    end

    log = NamedTuple{(:iter, :seconds, :objective, :residual, :status),
                     Tuple{Int,Float64,Float64,Float64,String}}[]
    result = nothing
    nvar = nbin = ncon = 0
    converged = false
    total_solve = 0.0

    for it in 1:max_iter
        model = Model(optimizer)
        silent && set_silent(model)
        for (k, val) in attributes
            set_optimizer_attribute(model, k, val)
        end
        time_limit_sec === nothing || set_time_limit_sec(model, time_limit_sec)

        # ---- host variables ---------------------------------------------------------
        @variable(model, 0 <= Pgen[SLACK_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Qgen[SLACK_SET, HOUR_SET, QUARTER_SET])
        @variable(model, c.Vmin <= v[BUS_SET, HOUR_SET, QUARTER_SET] <= c.Vmax, start = c.Vnom)
        @variable(model, v_r[BUS_SET, HOUR_SET, QUARTER_SET], start = c.Vnom)
        @variable(model, v_im[BUS_SET, HOUR_SET, QUARTER_SET], start = 0)
        @variable(model, Ibr_r[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Ibr_im[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Psnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Qsnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, 0 <= Ploss[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, 0 <= Qloss[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Ibs_r[BUS_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Ibs_im[BUS_SET, HOUR_SET, QUARTER_SET])
        @variable(model, 0 <= Pdg[DG_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Qdg[DG_SET, HOUR_SET, QUARTER_SET], start = 0)
        @variable(model, 0 <= PVC[DG_SET, HOUR_SET, QUARTER_SET])

        # ---- the droop module -------------------------------------------------------
        add_droop!(model, method, curve, v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar)

        # ---- inverter capability: a 16-segment polygon inscribing the S-circle ------
        k = 16
        for l in 1:k
            θ = l * π / k
            @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                cos(θ) * Pdg[d,h,m] + sin(θ) * Qdg[d,h,m] <=  c.Sdg_max[d])
            @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                cos(θ) * Pdg[d,h,m] + sin(θ) * Qdg[d,h,m] >= -c.Sdg_max[d])
        end
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Pdg[d,h,m] <= c.Pdg_max_vary[d][h,m])
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Pdg[d,h,m] <= c.Pdg_max[d])
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Qdg[d,h,m] <= c.Sdg_max[d])

        # ---- slack reference --------------------------------------------------------
        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET], v_r[i,h,m]  == c.Vnom)
        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET], v_im[i,h,m] == 0)

        # ---- Ohm's law along each branch (exact, linear) ----------------------------
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            v_r[i,h,m] - v_r[j,h,m] == R[(i,j)]*Ibr_r[(i,j),h,m] - X[(i,j)]*Ibr_im[(i,j),h,m])
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            v_im[i,h,m] - v_im[j,h,m] == R[(i,j)]*Ibr_im[(i,j),h,m] + X[(i,j)]*Ibr_r[(i,j),h,m])

        # ---- KCL at every bus (exact, linear) ---------------------------------------
        @constraint(model, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
            Ibs_r[bus,h,m] == sum(Ibr_r[(bus,j),h,m] for (i,j) in BRANCH_SET if i == bus)
                            - sum(Ibr_r[(i,bus),h,m] for (i,j) in BRANCH_SET if j == bus))
        @constraint(model, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
            Ibs_im[bus,h,m] == sum(Ibr_im[(bus,j),h,m] for (i,j) in BRANCH_SET if i == bus)
                             - sum(Ibr_im[(i,bus),h,m] for (i,j) in BRANCH_SET if j == bus))

        # ---- power balance: v·I linearised about the previous iterate ---------------
        Plin(i,h,m) = v_r_pr[i,h,m]*Ibs_r[i,h,m]  + v_im_pr[i,h,m]*Ibs_im[i,h,m] +
                      Ibs_r_pr[i,h,m]*v_r[i,h,m]  + Ibs_im_pr[i,h,m]*v_im[i,h,m] -
                      v_r_pr[i,h,m]*Ibs_r_pr[i,h,m] - v_im_pr[i,h,m]*Ibs_im_pr[i,h,m]
        Qlin(i,h,m) = v_im_pr[i,h,m]*Ibs_r[i,h,m] - v_r_pr[i,h,m]*Ibs_im[i,h,m] +
                      Ibs_r_pr[i,h,m]*v_im[i,h,m] - Ibs_im_pr[i,h,m]*v_r[i,h,m] -
                      v_im_pr[i,h,m]*Ibs_r_pr[i,h,m] + v_r_pr[i,h,m]*Ibs_im_pr[i,h,m]

        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Pgen[i,h,m] == Plin(i,h,m))
        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Qgen[i,h,m] == Qlin(i,h,m))
        @constraint(model, [i in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Pload[i,h,m] == Plin(i,h,m))
        @constraint(model, [i in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Qload[i,h,m] == Qlin(i,h,m))
        @constraint(model, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Pload[i,h,m] + Pdg[i,h,m] == Plin(i,h,m))
        @constraint(model, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Qload[i,h,m] + Qdg[i,h,m] == Qlin(i,h,m))

        # ---- branch loss: |I|² linearised about the previous iterate ----------------
        Isq(i,j,h,m) = 2*Ibr_r_pr[((i,j),h,m)]*Ibr_r[(i,j),h,m]   - Ibr_r_pr[((i,j),h,m)]^2 +
                       2*Ibr_im_pr[((i,j),h,m)]*Ibr_im[(i,j),h,m] - Ibr_im_pr[((i,j),h,m)]^2
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Ploss[(i,j),h,m] == R[(i,j)] * Isq(i,j,h,m))
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Qloss[(i,j),h,m] == X[(i,j)] * Isq(i,j,h,m))

        # ---- sending-end powers and branch losses ------------------------------------
        # Carrying the sending-end powers and the loss variables alongside the current
        # formulation is the extension to the base IVACOPF introduced in Emami Mirak &
        # Inaolaji (EPSR, doi:10.1016/j.epsr.2026.113613). It ties each bus injection to
        # the flows leaving it and each branch pair to its loss, which pins down the
        # dispatch more tightly than the current equations alone.
        @constraint(model, [bus in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Pgen[bus,h,m] == sum(Psnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Pload[bus,h,m])
        @constraint(model, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Pdg[bus,h,m] == sum(Psnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Pload[bus,h,m])
        @constraint(model, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            0 == sum(Psnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Pload[bus,h,m])
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Psnd[(i,j),h,m] == Ploss[(i,j),h,m] - Psnd[(j,i),h,m])
        @constraint(model, [bus in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Qgen[bus,h,m] == sum(Qsnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Qload[bus,h,m])
        @constraint(model, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Qdg[bus,h,m] == sum(Qsnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Qload[bus,h,m])
        @constraint(model, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            0 == sum(Qsnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Qload[bus,h,m])
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Qsnd[(i,j),h,m] == Qloss[(i,j),h,m] - Qsnd[(j,i),h,m])

        # ---- voltage magnitude, linearised about the previous iterate ---------------
        @constraint(model, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
            v[i,h,m] == (v_r_pr[i,h,m]/sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2))*v_r[i,h,m]
                      + (v_im_pr[i,h,m]/sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2))*v_im[i,h,m])

        # ---- curtailment and objective ----------------------------------------------
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            PVC[d,h,m] == c.Pdg_max_vary[d][h,m] - Pdg[d,h,m])
        @objective(model, Min, sum(PVC[d,h,m] for d in DG_SET, h in HOUR_SET, m in QUARTER_SET))

        nvar = num_variables(model)
        nbin = count(is_binary, all_variables(model))
        ncon = num_constraints(model; count_variable_in_set_constraints = false)

        secs = @elapsed optimize!(model)
        total_solve += secs
        status = termination_status(model)
        # ALMOST_LOCALLY_SOLVED is tolerated, and recorded: the Heaviside encoding is
        # non-smooth at every breakpoint, so an interior-point solver routinely stops
        # just short of its convergence tolerance. That is the price of the encoding,
        # not a failure of the model, and the tutorial reports it rather than hiding it.
        status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ||
            error("iteration $it terminated with status $status")

        # residual of the exact (bilinear) loss identity: how far the linearisation is
        # from the true power flow
        residual = maximum(abs(value((v_r[i,h,m] - v_r[j,h,m]) * Ibr_r[(i,j),h,m]
                                   + (v_im[i,h,m] - v_im[j,h,m]) * Ibr_im[(i,j),h,m]
                                   - Ploss[(i,j),h,m]))
                           for (i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET)

        push!(log, (iter = it, seconds = secs, objective = objective_value(model),
                    residual = residual, status = string(status)))

        # refresh the linearisation point
        for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET
            v_r_pr[b,h,q]    = value(v_r[b,h,q])
            v_im_pr[b,h,q]   = value(v_im[b,h,q])
            Ibs_r_pr[b,h,q]  = value(Ibs_r[b,h,q])
            Ibs_im_pr[b,h,q] = value(Ibs_im[b,h,q])
        end
        for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET
            Ibr_r_pr[(br,h,q)]  = value(Ibr_r[br,h,q])
            Ibr_im_pr[(br,h,q)] = value(Ibr_im[br,h,q])
        end

        result = (V   = [value(v[i,h,m])    for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
                  Vdg = [value(v[d,h,m])    for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  Qdg = [value(Qdg[d,h,m])  for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  Pdg = [value(Pdg[d,h,m])  for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  PVC = [value(PVC[d,h,m])  for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  Ploss = sum(value(Ploss[(i,j),h,m])
                              for (i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET))

        if residual < tol
            converged = true
            break
        end
    end

    # the warm-start solve is part of the cost of the answer, so it is counted
    return DOPFResult(result.V, result.Vdg, result.Qdg, result.Pdg, result.PVC,
                      result.Ploss, log, nvar, nbin, ncon, converged,
                      total_solve + warm_seconds)
end

# Expand an inverter dispatch into a full complex network state with an exact
# backward/forward sweep, and return it in the form the IVACOPF linearisation expects.
#
# The sweep solves the true AC power flow for the given injections, so the resulting
# point satisfies Ohm's law and KCL exactly, unlike the LinDistFlow solution it is
# derived from, which satisfies neither. That is the whole value of it as a starting
# point: `Pdg`/`Qdg` supply a sensible *dispatch*, and the sweep turns it into a
# physically consistent *state*.
function _sweep_state(c::Case, Pdg::Array{Float64,3}, Qdg::Array{Float64,3})
    nb, nbr = length(c.BUS_SET), length(c.BRANCH_SET)
    dgpos   = Dict(d => i for (i, d) in enumerate(c.DG_SET))

    v_r  = ones(nb, 24, 4);  v_im  = zeros(nb, 24, 4)
    Ibs_r = zeros(nb, 24, 4); Ibs_im = zeros(nb, 24, 4)
    Ibr_r  = Dict((br, h, q) => 0.0 for br in c.BRANCH_SET, h in c.HOUR_SET, q in c.QUARTER_SET)
    Ibr_im = Dict((br, h, q) => 0.0 for br in c.BRANCH_SET, h in c.HOUR_SET, q in c.QUARTER_SET)

    for h in c.HOUR_SET, q in c.QUARTER_SET
        # net consumption at each bus: load minus whatever the inverters inject
        Snet = map(c.BUS_SET) do b
            P = c.Pload[b,h,q]; Q = c.Qload[b,h,q]
            if haskey(dgpos, b)
                P -= Pdg[dgpos[b],h,q]; Q -= Qdg[dgpos[b],h,q]
            end
            ComplexF64(P, Q)
        end

        Vc  = ones(ComplexF64, nb)
        Ibr = zeros(ComplexF64, nbr)
        for _ in 1:40
            Ibus = [conj(Snet[b] / Vc[b]) for b in c.BUS_SET]
            for k in nbr:-1:1                                     # backward: currents
                (_, j) = c.BRANCH_SET[k]
                Ibr[k] = Ibus[j] + sum(Ibr[t] for t in 1:nbr
                                       if c.BRANCH_SET[t][1] == j; init = 0.0 + 0.0im)
            end
            for k in 1:nbr                                        # forward: voltages
                (i, j) = c.BRANCH_SET[k]
                Vc[j] = Vc[i] - (c.R[(i,j)] + im*c.X[(i,j)]) * Ibr[k]
            end
        end

        for b in c.BUS_SET
            v_r[b,h,q] = real(Vc[b]); v_im[b,h,q] = imag(Vc[b])
            # bus injection current, in the IVACOPF sign convention
            # (KCL there reads Ibs = outgoing − incoming, i.e. minus the load current)
            Ibus_b = conj(Snet[b] / Vc[b])
            Ibs_r[b,h,q] = -real(Ibus_b); Ibs_im[b,h,q] = -imag(Ibus_b)
        end
        for k in 1:nbr
            br = c.BRANCH_SET[k]
            Ibr_r[(br,h,q)] = real(Ibr[k]); Ibr_im[(br,h,q)] = imag(Ibr[k])
        end
    end
    return v_r, v_im, Ibs_r, Ibs_im, Ibr_r, Ibr_im
end

# The linearised branch-flow (LinDistFlow) host. Same droop module, same inverter
# constraints, same objective; only the network model differs, and there is no outer
# loop because nothing in it is nonlinear.
#
# Building and solving are split so the model can be inspected without a solver licence.
function _build_lindistflow(c::Case, optimizer;
                            method::Symbol, curve::DroopCurve, silent::Bool,
                            time_limit_sec, attributes)

    HOUR_SET, QUARTER_SET = c.HOUR_SET, c.QUARTER_SET
    BUS_SET, BRANCH_SET   = c.BUS_SET, c.BRANCH_SET
    DG_SET, NON_DG_SET    = c.DG_SET, c.NON_DG_SET
    SLACK_SET             = [c.slack]
    R, X                  = c.R, c.X
    qbar                  = Dict(d => c.Sdg_max[d] for d in DG_SET)

    model = Model(optimizer)
    silent && set_silent(model)
    for (k, val) in attributes
        set_optimizer_attribute(model, k, val)
    end
    time_limit_sec === nothing || set_time_limit_sec(model, time_limit_sec)

    @variable(model, Pgen[SLACK_SET, HOUR_SET, QUARTER_SET])
    @variable(model, Qgen[SLACK_SET, HOUR_SET, QUARTER_SET])
    @variable(model, c.Vmin <= v[BUS_SET, HOUR_SET, QUARTER_SET] <= c.Vmax, start = c.Vnom)
    @variable(model, Pbr[BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(model, Qbr[BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(model, 0 <= Pdg[DG_SET, HOUR_SET, QUARTER_SET])
    @variable(model, Qdg[DG_SET, HOUR_SET, QUARTER_SET], start = 0)
    @variable(model, 0 <= PVC[DG_SET, HOUR_SET, QUARTER_SET])

    # ---- the droop module: byte-for-byte the same call as the IVACOPF host ----------
    add_droop!(model, method, curve, v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar)

    # ---- inverter capability and limits (identical to the IVACOPF host) -------------
    k = 16
    for l in 1:k
        θ = l * π / k
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            cos(θ) * Pdg[d,h,m] + sin(θ) * Qdg[d,h,m] <=  c.Sdg_max[d])
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            cos(θ) * Pdg[d,h,m] + sin(θ) * Qdg[d,h,m] >= -c.Sdg_max[d])
    end
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Pdg[d,h,m] <= c.Pdg_max_vary[d][h,m])
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Pdg[d,h,m] <= c.Pdg_max[d])
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Qdg[d,h,m] <= c.Sdg_max[d])

    # ---- slack reference ------------------------------------------------------------
    @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET], v[i,h,m] == c.Vnom)

    # ---- power balance: net injection = outflow − inflow, losses neglected ----------
    outflow(P, b, h, m) = sum(P[(i,j),h,m] for (i,j) in BRANCH_SET if i == b; init = 0.0)
    inflow(P, b, h, m)  = sum(P[(i,j),h,m] for (i,j) in BRANCH_SET if j == b; init = 0.0)

    @constraint(model, [b in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
        Pgen[b,h,m] - c.Pload[b,h,m] == outflow(Pbr, b, h, m) - inflow(Pbr, b, h, m))
    @constraint(model, [b in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
        Qgen[b,h,m] - c.Qload[b,h,m] == outflow(Qbr, b, h, m) - inflow(Qbr, b, h, m))
    @constraint(model, [b in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
        -c.Pload[b,h,m] == outflow(Pbr, b, h, m) - inflow(Pbr, b, h, m))
    @constraint(model, [b in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
        -c.Qload[b,h,m] == outflow(Qbr, b, h, m) - inflow(Qbr, b, h, m))
    @constraint(model, [b in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Pdg[b,h,m] - c.Pload[b,h,m] == outflow(Pbr, b, h, m) - inflow(Pbr, b, h, m))
    @constraint(model, [b in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Qdg[b,h,m] - c.Qload[b,h,m] == outflow(Qbr, b, h, m) - inflow(Qbr, b, h, m))

    # ---- voltage drop along each branch ---------------------------------------------
    @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
        v[j,h,m] == v[i,h,m] -
            (R[(i,j)] * Pbr[(i,j),h,m] + X[(i,j)] * Qbr[(i,j),h,m]) / c.Vnom)

    # ---- curtailment and objective (identical to the IVACOPF host) ------------------
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        PVC[d,h,m] == c.Pdg_max_vary[d][h,m] - Pdg[d,h,m])
    @objective(model, Min, sum(PVC[d,h,m] for d in DG_SET, h in HOUR_SET, m in QUARTER_SET))

    return (; model, v, Pbr, Qbr, Pdg, Qdg, PVC)
end

function _solve_lindistflow(c::Case, optimizer; kwargs...)
    b = _build_lindistflow(c, optimizer; kwargs...)
    model, v, Pbr, Qbr = b.model, b.v, b.Pbr, b.Qbr
    Pdg, Qdg, PVC      = b.Pdg, b.Qdg, b.PVC
    BUS_SET, BRANCH_SET, DG_SET = c.BUS_SET, c.BRANCH_SET, c.DG_SET
    HOUR_SET, QUARTER_SET       = c.HOUR_SET, c.QUARTER_SET
    R = c.R

    nvar = num_variables(model)
    nbin = count(is_binary, all_variables(model))
    ncon = num_constraints(model; count_variable_in_set_constraints = false)

    secs = @elapsed optimize!(model)
    status = termination_status(model)
    status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ||
        error("LinDistFlow solve terminated with status $status")

    # Losses are not in the model. This is the standard post-hoc estimate from the
    # solved flows, reported so the field means something. It is not a model quantity.
    Ploss = sum(R[(i,j)] * (value(Pbr[(i,j),h,m])^2 + value(Qbr[(i,j),h,m])^2) /
                value(v[i,h,m])^2
                for (i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET)

    log = [(iter = 1, seconds = secs, objective = objective_value(model),
            residual = 0.0, status = string(status))]

    return DOPFResult(
        [value(v[i,h,m])   for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
        [value(v[d,h,m])   for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
        [value(Qdg[d,h,m]) for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
        [value(Pdg[d,h,m]) for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
        [value(PVC[d,h,m]) for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
        Ploss, log, nvar, nbin, ncon, true, secs)
end

"""
    base_case_voltages(case)

Voltage magnitudes for the same feeder and demand with **no** smart inverters, obtained
by a backward/forward sweep power flow. This is the reference the droop-aware dispatch
is compared against. The branch list of this feeder is already in topological order.
"""
function base_case_voltages(c::Case)
    nb, nbr = length(c.BUS_SET), length(c.BRANCH_SET)
    V = zeros(nb, 24, 4)
    for h in c.HOUR_SET, q in c.QUARTER_SET
        Vc  = ones(ComplexF64, nb)
        Ibr = zeros(ComplexF64, nbr)
        for _ in 1:30
            Ibus = [conj((c.Pload[b,h,q] + im*c.Qload[b,h,q]) / Vc[b]) for b in c.BUS_SET]
            for k in nbr:-1:1                                    # backward: currents
                (_, j) = c.BRANCH_SET[k]
                Ibr[k] = Ibus[j] + sum(Ibr[t] for t in 1:nbr
                                       if c.BRANCH_SET[t][1] == j; init = 0.0 + 0.0im)
            end
            for k in 1:nbr                                       # forward: voltages
                (i, j) = c.BRANCH_SET[k]
                Vc[j] = Vc[i] - (c.R[(i,j)] + im*c.X[(i,j)]) * Ibr[k]
            end
        end
        V[:,h,q] = abs.(Vc)
    end
    return V
end
