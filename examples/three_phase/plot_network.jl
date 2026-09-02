# Schematic of network_5_Feeder_2.
#
# The feeder is 194 buses but only ~40 of them are interesting: the rest are degree-2
# points in the middle of a cable run. This collapses those chains, so what is drawn is
# the source, every branch point, every load, every PV site and every dead end, joined
# by edges whose horizontal length is the real electrical distance.
#
#   julia --project=. plot_network.jl

using JSON3, Printf, Plots

const DATA = joinpath(@__DIR__, "data")
const CASE = "network_5_Feeder_2"
const N_PV_PER_PHASE = 4          # keep in step with LinDist3Flow_Lambda.jl

net   = JSON3.read(read(joinpath(DATA, "$CASE.bmopf.json"), String))
BUSES = [String(k) for k in keys(net.bus)]
bid   = Dict(b => i for (i, b) in enumerate(BUSES));  nb = length(BUSES)
SLACK = String(net.voltage_source.source.bus)
VBASE = Float64(net.voltage_source.source.v_magnitude[1])
SBASE = 100e3;  ZBASE = VBASE^2 / SBASE

# ---- topology and per-branch impedance magnitude --------------------------------------
F = String[]; T = String[]; ZM = Float64[]
for (_, l) in pairs(net.line)
    lc = net.linecode[Symbol(l.linecode)]
    z  = complex(Float64(lc[Symbol("R_series_1_1")]), Float64(lc[Symbol("X_series_1_1")]))
    push!(F, String(l.bus_from)); push!(T, String(l.bus_to))
    push!(ZM, abs(z) * Float64(l.length))          # ohms
end

adj = Dict(b => Tuple{String,Int}[] for b in BUSES)
for k in eachindex(F)
    push!(adj[F[k]], (T[k], k)); push!(adj[T[k]], (F[k], k))
end

# rooted tree: distance in ohms, parent, children, BFS order
DIST = Dict(SLACK => 0.0); PAR = Dict{String,String}()
KIDS = Dict(b => String[] for b in BUSES); ORDER = [SLACK]; q = [SLACK]
while !isempty(q)
    b = popfirst!(q)
    for (nbb, k) in adj[b]
        if !haskey(DIST, nbb)
            DIST[nbb] = DIST[b] + ZM[k]; PAR[nbb] = b
            push!(KIDS[b], nbb); push!(ORDER, nbb); push!(q, nbb)
        end
    end
end

# ---- loads and PV sites ---------------------------------------------------------------
loadkW = Dict{String,Float64}(); loadph = Dict{String,Int}()
for (_, ld) in pairs(net.load)
    b = String(ld.bus); φ = parse(Int, String(ld.terminal_map[1]))
    loadkW[b] = get(loadkW, b, 0.0) + Float64(ld.p_nom[1]) / 1e3
    loadph[b] = φ
end
# PV sites and their class: the same placement and class rotation as the OPF script
const PV_CLASSES = (("A", 3.0), ("B", 5.0), ("C", 8.0), ("D", 12.0))
PVPH  = Dict{String,Int}()      # bus → phase
PVCLS = Dict{String,Int}()      # bus → class index
for φ in 1:3
    cand = [b for b in BUSES if get(loadph, b, 0) == φ]
    sort!(cand, by = b -> -DIST[b])
    for (rank, b) in enumerate(cand[1:min(N_PV_PER_PHASE, length(cand))])
        PVPH[b]  = φ
        PVCLS[b] = mod(rank - 1 + (φ - 1), length(PV_CLASSES)) + 1
    end
end
PVBUS = PVPH   # kept for the `iskey` test below

# ---- keep only the structurally interesting buses, contract the rest -------------------
iskey(b) = b == SLACK || haskey(loadkW, b) || haskey(PVBUS, b) ||
           length(KIDS[b]) != 1
KEY = [b for b in ORDER if iskey(b)]
keyset = Set(KEY)

# each key bus links back to its nearest key ancestor
kpar = Dict{String,String}()
for b in KEY
    b == SLACK && continue
    p = PAR[b]
    while !(p in keyset); p = PAR[p]; end
    kpar[b] = p
end
kkids = Dict(b => String[] for b in KEY)
for b in KEY; b == SLACK || push!(kkids[kpar[b]], b); end

# ---- dendrogram layout: x = electrical distance, y = leaf ordering ---------------------
Y = Dict{String,Float64}(); leaf = Ref(0.0)
function place(b)
    if isempty(kkids[b])
        Y[b] = leaf[]; leaf[] += 1.0
    else
        for c in kkids[b]; place(c); end
        Y[b] = sum(Y[c] for c in kkids[b]) / length(kkids[b])
    end
end
place(SLACK)

@printf("full feeder : %d buses, %d lines\n", nb, length(F))
@printf("schematic   : %d nodes after collapsing series chains\n", length(KEY))
@printf("branch pts  : %d | loads: %d | PV sites: %d | dead ends: %d\n",
        count(b -> length(KIDS[b]) > 1, BUSES), length(loadkW), length(PVBUS),
        count(b -> isempty(KIDS[b]), BUSES))

# ---- draw ------------------------------------------------------------------------------
# Colour carries the inverter class, matching the droop figure. Phase is written into
# each label, so the two do not fight for the same visual channel.
CLSCOL = [:seagreen, :orangered, :dodgerblue, :mediumorchid]
CLSMS  = [9, 11, 13, 15]                       # marker size tracks the rating

xmax = maximum(values(DIST))
p = plot(size = (1280, 760), legend = :topright, framestyle = :box,
         grid = true, gridalpha = 0.12, legendfontsize = 9,
         xlabel = "electrical distance from the substation  (Ω)", ylabel = "",
         yticks = false,
         title = "network_5_Feeder_2 — 194-bus radial LV feeder, 12 smart inverters",
         titlefontsize = 12,
         xlims = (-0.008, xmax * 1.30),         # headroom for the bus labels on the right
         ylims = (-2.5, maximum(values(Y)) + 9),  # headroom so the legend clears the tree
         left_margin = 5Plots.mm, bottom_margin = 6Plots.mm,
         right_margin = 5Plots.mm, top_margin = 6Plots.mm)

for b in KEY                                   # edges, elbow-style
    b == SLACK && continue
    a = kpar[b]
    plot!(p, [DIST[a], DIST[a], DIST[b]], [Y[a], Y[b], Y[b]],
          lw = 1.6, color = :grey60, label = false)
end

plain = [b for b in KEY if !haskey(loadkW, b) && b != SLACK]
scatter!(p, [DIST[b] for b in plain], [Y[b] for b in plain],
         m = :circle, ms = 2.4, mc = :grey55, msc = :grey55, label = "junction / dead end")

bare = [b for b in KEY if haskey(loadkW, b) && !haskey(PVPH, b)]   # load, no inverter
scatter!(p, [DIST[b] for b in bare], [Y[b] for b in bare], m = :circle, ms = 5.5,
         mc = :grey25, msc = :white, label = "load, no inverter")

for ci in eachindex(PV_CLASSES)                # PV + SI, coloured and sized by class
    bs = [b for b in KEY if get(PVCLS, b, 0) == ci]
    isempty(bs) && continue
    scatter!(p, [DIST[b] for b in bs], [Y[b] for b in bs],
             m = :star5, ms = CLSMS[ci], mc = CLSCOL[ci], msc = :black,
             label = "PV + SI, class $(PV_CLASSES[ci][1]), $(PV_CLASSES[ci][2]) kW")
end
for (b, ci) in PVCLS                           # bus number and phase beside each star
    annotate!(p, DIST[b] + 0.005, Y[b],
              text("$b  $(PV_CLASSES[ci][1])/φ$(PVPH[b])", 7, :left, :black))
end

scatter!(p, [DIST[SLACK]], [Y[SLACK]], m = :square, ms = 11, mc = :black,
         label = "substation  (240 V L-N)")
annotate!(p, xmax * 1.28, 1.5,
          text("vertical position is drawing order only; horizontal position is real",
               7, :right, :grey40))

savefig(p, joinpath(@__DIR__, "network_schematic.png"))
display(p)
println("\nwrote network_schematic.png")
