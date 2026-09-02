# Scalability of the three droop encodings on LinDist3Flow.
#
#   julia --project=. scalability.jl
#
# Runs each encoding on a small and a large real feeder, at the full 96-step horizon, and
# writes docs/src/assets/results/threephase/scalability.json.
#
# Nothing about the model changes between runs: the same three scripts are executed with
# TP_CASE / TP_STEPS / TP_NPV set in the environment. The question being asked is whether
# the droop encodings survive a twenty-fold larger network, not whether a bespoke
# formulation can be made to.

using JSON3, Printf

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "..", ".."))
const OUT  = joinpath(ROOT, "docs", "src", "assets", "results", "threephase")
mkpath(OUT)

const CASES = [(case = "network_5_Feeder_2",  label = "194-bus",  npv = 4),
               (case = "network_17_Feeder_6", label = "3856-bus", npv = 4)]
const SCRIPTS = [("lambda", "LinDist3Flow_Lambda.jl"),
                 ("bigm",   "LinDist3Flow_BigM.jl"),
                 ("heaviside", "LinDist3Flow_Heaviside.jl")]
const TIMEOUT = parse(Int, get(ENV, "TP_TIMEOUT", "900"))   # seconds per run
# horizons to fall back to, in order, when the full 96-step run does not finish
const REDUCED_STEPS = [parse(Int, x) for x in split(get(ENV, "TP_REDUCED", "12"), ',')]

"Run one script in a fresh process and scrape the numbers it prints."
function run_one(script, case, npv, steps)
    env = copy(ENV)
    env["TP_CASE"] = case; env["TP_NPV"] = string(npv); env["TP_STEPS"] = string(steps)
    env["GKSwstype"] = "100"
    cmd = Cmd(`$(Base.julia_cmd()) --project=$HERE $(joinpath(HERE, script))`, env = env)
    tmp = tempname()
    t0  = time()
    # A run that will not finish is a scalability result too, so bound it rather than
    # letting one configuration stall the whole sweep.
    proc = run(pipeline(cmd, stdout = tmp, stderr = devnull), wait = false)
    timedout = false
    while process_running(proc)
        if time() - t0 > TIMEOUT
            timedout = true; kill(proc); break
        end
        sleep(2)
    end
    wait(proc)
    wall = time() - t0
    ok  = !timedout && success(proc)
    txt = isfile(tmp) ? read(tmp, String) : ""
    rm(tmp, force = true)
    grab(re, default = NaN) = (m = match(re, txt)) === nothing ? default : parse(Float64, m.captures[1])
    iv(x) = isfinite(x) ? Int(x) : 0        # a killed run prints nothing to scrape
    return (; ok, timedout, wall,
            nvar  = iv(grab(r"model\s+:\s+(\d+) variables")),
            nbin  = iv(grab(r"variables \((\d+) binary\)")),
            ncon  = iv(grab(r"binary\), (\d+) constraints")),
            solve = grab(r"solve time\s+:\s+([\d.]+) s"),
            curt  = grab(r"PV curtailment\s+:\s+([\d.]+) kWh"),
            dev   = grab(r"q_dispatch − q_curve\|\s*:\s*([\d.e\-+]+)"),
            vlo   = grab(r"voltage range\s+:\s+([\d.]+) –"),
            vhi   = grab(r"voltage range\s+:\s+[\d.]+ – ([\d.]+)"))
end

rows = []
const LOG = joinpath(OUT, "scalability_log.txt")
logline(s) = (open(LOG, "a") do io; println(io, s); end)

fin(x) = isfinite(x) ? x : nothing
report(tag, c, steps, r, note) = begin
    row = Dict("encoding" => tag, "feeder" => c.label, "case" => c.case,
               "n_si" => 3 * c.npv, "steps" => steps, "ok" => r.ok,
               "timed_out" => r.timedout,
               "nvar" => r.nvar, "nbin" => r.nbin, "ncon" => r.ncon,
               "solve_seconds" => fin(r.solve), "wall_seconds" => fin(r.wall),
               "curt_kWh" => fin(r.curt), "max_droop_deviation" => fin(r.dev),
               "Vmin" => fin(r.vlo), "Vmax" => fin(r.vhi))
    note === nothing || (row["note"] = note)
    push!(rows, row)
    line = @sprintf("%-11s %-9s %-4d %-6d %-10s %10d %8d %10.1f %9.1f %10.1e",
        tag, c.label, 3 * c.npv, steps,
        r.ok ? "solved" : (r.timedout ? "timeout" : "FAILED"),
        r.nvar, r.nbin, r.solve, r.wall, r.dev)
    println(line); flush(stdout); logline(line)
end

@printf("%-11s %-9s %-4s %-6s %-10s %10s %8s %10s %9s %10s
",
        "encoding", "feeder", "SIs", "steps", "status", "variables", "binaries",
        "solve (s)", "wall (s)", "max|Δq|")
for c in CASES, (tag, script) in SCRIPTS
    r = run_one(script, c.case, c.npv, 96)
    report(tag, c, 96, r, nothing)

    # A run that does not finish at the full horizon is a scalability result, not a dead
    # end: shrink the day and ask again. That is how the reduced-horizon Heaviside row on
    # the large feeder is produced, and it is what shows the encoding still *solves*,
    # reproducing the droop to round-off, once the model is small enough to differentiate.
    if !r.ok
        for steps in REDUCED_STEPS
            r2 = run_one(script, c.case, c.npv, steps)
            report(tag, c, steps, r2,
                   "reduced horizon; the 96-step model does not finish")
            r2.ok && break
        end
    end
end

open(joinpath(OUT, "scalability.json"), "w") do io
    JSON3.pretty(io, Dict("runs" => rows, "horizon" => 96))
end
println("\nwrote ", joinpath(OUT, "scalability.json"))
