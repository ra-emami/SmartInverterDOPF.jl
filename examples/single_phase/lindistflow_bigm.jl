# =====================================================================================
#  Volt-VAr droop in a distribution OPF
#  HOST   : LinDistFlow, a linearised branch-flow model, solved once
#  DROOP  : Big-M (MILP, solved with Gurobi)
#
#  IEEE 33-bus feeder, 24 h at 15-minute resolution, three PV smart inverters at buses
#  7, 18 and 33, minimising total PV curtailment.
#
#  This is a thin wrapper: the model lives in the package, and `method` and `host` are
#  independent choices. Every one of the six combinations is a script like this one.
#
#  Run:  julia --project=examples/single_phase examples/single_phase/lindistflow_bigm.jl
# =====================================================================================

using SmartInverterDOPF, Gurobi, Printf

case = load_case()

res = solve_dopf(case, Gurobi.Optimizer;
                 method     = :bigm,
                 host       = :lindistflow,
                 attributes = Dict{String,Any}("MIPGap" => 1e-3, "OutputFlag" => 0))

# ---- results ---------------------------------------------------------------------------
E_avail = kWh(case, sum(case.Pdg_max_vary[d][h, m]
                        for d in case.DG_SET, h in 1:24, m in 1:4))
E_curt  = kWh(case, sum(res.PVC))

println("\n================ LinDistFlow + Big-M ================")
@printf("model               : %d variables (%d binary), %d constraints\n",
        res.nvar, res.nbin, res.ncon)
@printf("passes              : %d %s\n", length(res.iterations),
        res.converged ? "(converged)" : "(did NOT converge)")
@printf("solve time          : %.1f s\n", res.solve_seconds)
@printf("PV energy available : %.2f kWh\n", E_avail)
@printf("PV curtailment      : %.2f kWh  (%.3f %%)\n", E_curt, 100 * E_curt / E_avail)
@printf("network loss        : %.2f kWh\n", kWh(case, res.Ploss))
@printf("voltage range       : %.4f - %.4f p.u.\n", minimum(res.V), maximum(res.V))

# ---- does the dispatch lie on the droop curve? ------------------------------------------
curve = ieee1547_curve()
dev = maximum(abs(res.Qdg[i, h, m] -
                  droop_q(curve, res.Vdg[i, h, m], case.Sdg_max[case.DG_SET[i]]))
              for i in eachindex(case.DG_SET), h in 1:24, m in 1:4)
@printf("max |q_dispatch - q_curve| : %.3e p.u.\n", dev)
