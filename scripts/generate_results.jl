# Regenerate the precomputed results the documentation is built from.
#
#   julia --project=scripts scripts/generate_results.jl           # Gurobi + Ipopt
#   julia --project=scripts scripts/generate_results.jl highs     # HiGHS  + Ipopt
#   julia --project=scripts scripts/generate_results.jl glpk      # GLPK   + Ipopt
#
# Writes one JSON per method into docs/src/assets/results/. The documentation reads
# those files and redraws every figure and table at build time, so building the docs
# needs no optimisation solver at all, only Plots and JSON3.
#
# NOTE: use Gurobi. It is free for academic users (see the tutorial's solver note) and
# fast on this model, and the committed results are generated with it. The HiGHS and GLPK
# options are kept only for experimentation: on this case one returned an infeasible
# status inside the successive-linearisation loop and the other was too slow to finish,
# so neither produces usable results here.

using SmartInverterDOPF
using JSON3

const ROOT   = dirname(@__DIR__)
const OUTDIR = joinpath(ROOT, "docs", "src", "assets", "results")
mkpath(OUTDIR)

milp = lowercase(get(ARGS, 1, "gurobi"))
if milp == "highs"
    using HiGHS
    milp_opt, milp_name = HiGHS.Optimizer, "HiGHS"
    milp_attrs = Dict{String,Any}("mip_rel_gap" => 0.001, "output_flag" => false)
    @warn "HiGHS is not expected to complete the linearisation loop on this case; see the header comment."
elseif milp == "glpk"
    using GLPK
    milp_opt, milp_name = GLPK.Optimizer, "GLPK"
    milp_attrs = Dict{String,Any}("mip_gap" => 0.001, "msg_lev" => 0)
    @warn "GLPK is not expected to finish even the first pass on this case; see the header comment."
else
    using Gurobi
    milp_opt, milp_name = Gurobi.Optimizer, "Gurobi"
    milp_attrs = Dict{String,Any}("MIPGap" => 0.001, "OutputFlag" => 0)
end
using Ipopt

case  = load_case(joinpath(ROOT, "data"))
curve = ieee1547_curve()
qbar  = Dict(d => case.Sdg_max[d] for d in case.DG_SET)

# time labels, straight from the profile file
tlabels = JSON3.read(read(joinpath(ROOT, "data", "load_profiles_15min.json"), String)).time

# the no-inverter reference, identical for all three methods
@info "solving the base case (no smart inverters)"
Vbase_prof = base_case_voltages(case)

# irradiance ceiling per inverter, flattened to 96 steps, in kW
avail_kW(d) = vec(permutedims(case.Pdg_max_vary[d])) .* (case.Sbase / 1e3)

methods = [(:bigm,      "bigm",      milp_opt,      milp_name, milp_attrs),
           (:lambda,    "lambda",    milp_opt,      milp_name, milp_attrs),
           (:heaviside, "heaviside", Ipopt.Optimizer, "Ipopt",
            Dict{String,Any}("max_iter" => 3000, "print_level" => 0))]

summary = Dict{String,Any}()

for (method, slug, opt, solver, attrs) in methods
    @info "solving" method solver
    res = solve_dopf(case, opt; method = method, curve = curve, attributes = attrs)

    E_avail = kWh(case, sum(sum(case.Pdg_max_vary[d]) for d in case.DG_SET))
    E_curt  = kWh(case, sum(res.PVC))

    # how far each dispatch point strays from the exact droop law
    dev = maximum(abs(res.Qdg[i,h,m] - droop_q(curve, res.Vdg[i,h,m], qbar[case.DG_SET[i]]))
                  for i in 1:ndg(case), h in case.HOUR_SET, m in case.QUARTER_SET)

    # flatten [dg, hour, quarter] into a 96-step series per inverter
    series(A, i, scale) = vec(permutedims(A[i, :, :])) .* scale
    kW = case.Sbase / 1e3

    payload = Dict{String,Any}(
        "method"        => slug,
        "solver"        => solver,
        "converged"     => res.converged,
        "solve_seconds" => res.solve_seconds,
        "n_iterations"  => length(res.iterations),
        "iterations"    => [Dict("iter" => r.iter, "seconds" => r.seconds,
                                 "objective" => r.objective, "residual" => r.residual,
                                 "status" => r.status) for r in res.iterations],
        "nvar" => res.nvar, "nbin" => res.nbin, "ncon" => res.ncon,
        "E_avail_kWh"  => E_avail,
        "E_curt_kWh"   => E_curt,
        "E_deliv_kWh"  => E_avail - E_curt,
        "curt_percent" => 100 * E_curt / E_avail,
        "loss_kWh"     => kWh(case, res.Ploss),
        "Vmin" => minimum(res.V), "Vmax" => maximum(res.V),
        "max_droop_deviation" => dev,
        # scatter of every inverter operating point against the curve
        "Vdg_flat"     => vec(res.Vdg),
        "Qdg_norm"     => vec(res.Qdg ./ [qbar[d] for d in case.DG_SET,
                                          h in case.HOUR_SET, m in case.QUARTER_SET]),
        # per-bus voltage envelope over the whole day
        "V_with_max" => vec(maximum(res.V, dims = (2,3))),
        "V_with_min" => vec(minimum(res.V, dims = (2,3))),
        # per-inverter daily series, 96 steps
        "Vdg_series" => [series(res.Vdg, i, 1.0)  for i in 1:ndg(case)],
        "Qdg_series" => [series(res.Qdg, i, kW)   for i in 1:ndg(case)],
        "Pdg_series" => [series(res.Pdg, i, kW)   for i in 1:ndg(case)],
        "PVC_series" => [series(res.PVC, i, kW)   for i in 1:ndg(case)],
    )
    open(joinpath(OUTDIR, "$slug.json"), "w") do io
        JSON3.pretty(io, payload)
    end
    summary[slug] = payload
    @info "  done" curtail_kWh=E_curt iterations=length(res.iterations) deviation=dev
end

# everything shared by the three method pages
meta = Dict{String,Any}(
    "n_bus"        => nbus(case),
    "n_branch"     => length(case.BRANCH_SET),
    "n_steps"      => 96,
    "Sbase_kVA"    => case.Sbase / 1e3,
    "Vbase_kV"     => case.Vbase / 1e3,
    "time"         => collect(tlabels),
    "DG_SET"       => case.DG_SET,
    "Pdg_max_kW"   => [case.Pdg_max[d] * case.Sbase / 1e3 for d in case.DG_SET],
    "Sdg_max_kVA"  => [case.Sdg_max[d] * case.Sbase / 1e3 for d in case.DG_SET],
    "qbar_kVAr"    => [qbar[d] * case.Sbase / 1e3 for d in case.DG_SET],
    "avail_kW"     => [avail_kW(d) for d in case.DG_SET],
    "Vbp"          => curve.Vbp,
    "qshape"       => curve.qshape,
    "Vmin_limit"   => case.Vmin,
    "Vmax_limit"   => case.Vmax,
    "V_base_max"   => vec(maximum(Vbase_prof, dims = (2,3))),
    "V_base_min"   => vec(minimum(Vbase_prof, dims = (2,3))),
    "V_base_lo"    => minimum(Vbase_prof),
    "V_base_hi"    => maximum(Vbase_prof),
    "julia_version" => string(VERSION),
)
open(joinpath(OUTDIR, "case.json"), "w") do io
    JSON3.pretty(io, meta)
end

println("\nwrote $(OUTDIR)")
for slug in ("bigm", "lambda", "heaviside")
    p = summary[slug]
    println(rpad(slug, 11), " ", rpad(p["solver"], 7),
            " iters ", p["n_iterations"],
            " | curtail ", round(p["E_curt_kWh"], digits = 3), " kWh",
            " | loss ", round(p["loss_kWh"], digits = 2), " kWh",
            " | V ", round(p["Vmin"], digits = 4), "-", round(p["Vmax"], digits = 4))
end
