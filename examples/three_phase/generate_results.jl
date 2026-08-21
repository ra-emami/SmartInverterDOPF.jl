# Regenerate the precomputed three-phase results the documentation is built from.
#
#   julia --project=. examples/three_phase/generate_results.jl
#
# Each of the six shipped scripts is executed in its own module and its results are
# harvested, so the numbers in the documentation come from exactly the code a reader
# runs — there is no second copy of the model here to drift out of step.
#
# Writes docs/src/assets/results/threephase/
#   case.json                                     — feeder, fleet and droop description
#   {lambda,bigm,heaviside}.json                  — LinDist3Flow host
#   iva_{lambda,bigm,heaviside}.json              — IVACOPF host
#
# Set TP_HOSTS=lindist3flow or TP_HOSTS=ivacopf to regenerate only one family.
# case.json is rewritten only when absent — delete it to refresh the case description.

using JSON3, Printf

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "..", ".."))
const OUT  = joinpath(ROOT, "docs", "src", "assets", "results", "threephase")
mkpath(OUT)

const HOSTS = split(get(ENV, "TP_HOSTS", "lindist3flow,ivacopf"), ',')

"Run one of the standalone scripts in a private module and hand back its bindings."
function run_script(file)
    @info "running" file
    m = Module(Symbol("TP_", file))
    Core.eval(m, :(using Base))
    Base.include(m, joinpath(HERE, file))
    return m
end

# ---- the exact AC audit, applied to a solved module ------------------------------------
#  Both hosts expose the same names, so one audit serves both: `sweep(t)` is an exact
#  three-phase backward/forward power flow on the solved dispatch, and the question it
#  answers is whether the voltage the host told each inverter to read is the voltage it
#  would really see. `dev_model` is the exactness of the droop *encoding*; `dev_true` is
#  the accuracy of the *host*, and the two are entirely separate things.
function audit(m)
    dq(v, qb) = v <= m.VBP[2] ? qb :
                v <= m.VBP[3] ? qb * (m.VBP[3] - v) / (m.VBP[3] - m.VBP[2]) :
                v <= m.VBP[4] ? 0.0 :
                v <= m.VBP[5] ? -qb * (v - m.VBP[4]) / (m.VBP[5] - m.VBP[4]) : -qb
    gap = 0.0; dev_true = 0.0; dev_model = 0.0; nviol = 0
    tmin = Inf; tmax = -Inf
    for t in 1:m.T
        Vt = m.sweep(t)
        gap  = max(gap, maximum(abs.(m.V[:, :, t] .- Vt)))
        tmin = min(tmin, minimum(Vt)); tmax = max(tmax, maximum(Vt))
        nviol += count(<(m.VLIM[1] - 1e-6), Vt) + count(>(m.VLIM[2] + 1e-6), Vt)
        for i in 1:m.npv
            b, φ, qb = m.PV[i].bus, m.PV[i].phase, m.PV[i].Smax
            dev_true  = max(dev_true,  abs(m.Qdg_v[i, t] - dq(Vt[b, φ], qb)))
            dev_model = max(dev_model, abs(m.Qdg_v[i, t] - dq(m.V[b, φ, t], qb)))
        end
    end
    return (; gap, dev_true, dev_model, nviol, true_lo = tmin, true_hi = tmax)
end

# ---- the payload shared by both hosts ---------------------------------------------------
function base_payload(m, tag, host)
    kW(x)   = x * m.SBASE / 1e3
    E_avail = m.kWh(sum(m.Pavail)); E_curt = m.kWh(sum(m.PVC_v))
    a       = audit(m)
    return Dict(
        "host" => host,
        "method" => m.METHOD,
        "solver" => tag == "heaviside" ? "Ipopt" : "Gurobi",
        "model_class" => tag == "heaviside" ? "NLP" : "MILP",
        "nvar" => m.num_variables(m.model),
        "nbin" => count(m.is_binary, m.all_variables(m.model)),
        "ncon" => m.num_constraints(m.model; count_variable_in_set_constraints = false),
        "status" => string(m.status),
        "E_avail_kWh" => E_avail, "E_curt_kWh" => E_curt,
        "curt_percent" => 100 * E_curt / E_avail,
        "Vmin" => minimum(m.V), "Vmax" => maximum(m.V),
        "Vmin_phase" => [minimum(m.V[:, φ, :]) for φ in 1:3],
        "Vmax_phase" => [maximum(m.V[:, φ, :]) for φ in 1:3],
        # per-time-step envelopes, for the daily figure
        "Vmax_t" => [[maximum(m.V[:, φ, t]) for t in 1:m.T] for φ in 1:3],
        "Vmin_t" => [[minimum(m.V[:, φ, t]) for t in 1:m.T] for φ in 1:3],
        # per-inverter series, for the droop-verification figure
        "Vdg_series" => [[m.V[m.PV[i].bus, m.PV[i].phase, t] for t in 1:m.T] for i in 1:m.npv],
        "Qdg_series" => [[m.Qdg_v[i, t] for t in 1:m.T] for i in 1:m.npv],
        "P_avail_kW" => [kW(sum(m.Pavail[:, t])) for t in 1:m.T],
        "P_disp_kW"  => [kW(sum(m.Pdg_v[:, t]))  for t in 1:m.T],
        # verification
        "max_droop_deviation" => a.dev_model,
        "audit" => Dict("v_gap" => a.gap, "droop_residual_true_v" => a.dev_true,
                        "n_limit_violations" => a.nviol,
                        "true_Vmin" => a.true_lo, "true_Vmax" => a.true_hi),
    )
end

# ---- case description, taken from whichever module ran first ----------------------------
function write_case(L)
    kW(x) = x * L.SBASE / 1e3
    case = Dict(
        "name" => L.CASE,
        "n_bus" => L.nb, "n_line" => L.nbr, "n_load" => length(L.net.load),
        "Vbase_V" => L.VBASE, "Sbase_kVA" => L.SBASE_KVA, "Zbase_ohm" => L.ZBASE,
        "length_m" => sum(Float64(l.length) for (_, l) in pairs(L.net.line)),
        "loads_per_phase" => [count(b -> L.Pload_pk[b, φ] > 0, 1:L.nb) for φ in 1:3],
        "load_kW_per_phase" => [kW(sum(L.Pload_pk[:, φ])) for φ in 1:3],
        "load_kW_total" => kW(sum(L.Pload_pk)),
        "Vmin_limit" => L.VLIM[1], "Vmax_limit" => L.VLIM[2],
        "Vbp" => L.VBP, "qshape" => L.QSHAPE,
        "n_steps" => L.T,
        "classes" => [Dict("name" => c[1], "P_kW" => c[2],
                           "S_kVA" => L.S_OVER_P * c[2],
                           "qbar_pu" => L.S_OVER_P * c[2] * 1e3 / L.SBASE)
                      for c in L.PV_CLASSES],
        "sites" => [Dict("bus" => L.BUSES[g.bus], "phase" => g.phase,
                         "class" => L.PV_CLASSES[g.cls][1],
                         "class_idx" => g.cls,
                         "Z_ohm" => L.DIST[L.BUSES[g.bus]] * L.ZBASE,
                         "P_kW" => kW(g.Pmax), "S_kVA" => kW(g.Smax))
                    for g in L.PV],
        "PV_kW_total" => sum(kW(g.Pmax) for g in L.PV),
    )
    open(joinpath(OUT, "case.json"), "w") do io; JSON3.pretty(io, case); end
end

const TAGS = ("lambda", "bigm", "heaviside")
const SUFFIX = Dict("lambda" => "Lambda", "bigm" => "BigM", "heaviside" => "Heaviside")

# ---- emit one solved module -------------------------------------------------------------
#  Everything that reads a freshly-included module lives in here, and the callers reach it
#  through `Base.invokelatest`. Under Julia 1.12's stricter binding world-age rules, code
#  compiled before `Base.include` ran cannot see the globals that `include` created; the
#  symptom is `UndefVarError: nb not defined ... binding may be too new`. Running the
#  accessor in the latest world is the documented fix.
function emit(m, tag, host)
    p = base_payload(m, tag, host)
    if host == "IVACOPF"
        # the successive-linearisation loop is what is extra here
        p["solve_seconds"]     = m.total_solve
        p["last_pass_seconds"] = m.secs
        p["warm_seconds"]      = m.warm_seconds
        p["warm_start"]        = m.WARMSTART ? "sweep" : "flat"
        p["n_passes"]          = length(m.iterlog)
        p["converged"]         = m.converged
        p["tol"]               = m.TOL
        p["loss_kWh"]          = m.loss_kWh
        p["peak_line_loading"] = m.Iutil
        p["iterations"] = [Dict("iter" => r.iter, "seconds" => r.seconds,
                                "objective" => r.objective, "MAPB" => r.MAPB,
                                "MRPB" => r.MRPB, "MVM" => r.MVM, "status" => r.status)
                           for r in m.iterlog]
    else
        p["solve_seconds"] = m.secs
    end
    name = host == "IVACOPF" ? "iva_$tag" : tag
    open(joinpath(OUT, "$name.json"), "w") do io; JSON3.pretty(io, p); end
    @printf("%-12s %-10s %-4s %6.1f s / %d pass%s  curtail %6.2f kWh (%.3f %%)  V %.4f\u2013%.4f  dev %.2e  true-v dev %.2e\n",
            host, tag, p["model_class"], p["solve_seconds"], get(p, "n_passes", 1),
            get(p, "n_passes", 1) == 1 ? " " : "es", p["E_curt_kWh"], p["curt_percent"],
            p["Vmin"], p["Vmax"], p["max_droop_deviation"],
            p["audit"]["droop_residual_true_v"])
    return p
end

# ---- run one host family ------------------------------------------------------------------
# `run_script` includes each script, so every module access afterwards goes through
# `Base.invokelatest` — see the note on `emit`.
function run_family(prefix, host)
    mods = Dict{String,Module}()
    for t in TAGS
        mods[t] = run_script(prefix * SUFFIX[t] * ".jl")
    end
    isfile(joinpath(OUT, "case.json")) ||
        Base.invokelatest(write_case, mods["lambda"])
    for t in TAGS
        Base.invokelatest(emit, mods[t], t, host)
    end
    return mods
end

"lindist3flow" in HOSTS && run_family("LinDist3Flow_", "LinDist3Flow")
"ivacopf"      in HOSTS && run_family("IVACOPF3Ph_",   "IVACOPF")

println("\nwrote ", OUT)
