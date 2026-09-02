# The IEEE 1547 Volt-VAr (Q-V) droop, and the three ways of writing it as
# constraints a solver will accept.

"""
    DroopCurve(Vbp, qshape)

A five-segment piecewise-linear Volt-VAr characteristic, given by six breakpoint
voltages `Vbp` (increasing, p.u.) and the reactive set-point at each breakpoint
`qshape`, expressed as a multiple of the inverter's reactive capability `q̄`.

The IEEE 1547 default shape is `[1, 1, 0, 0, -1, -1]`: full injection below the
lower knee, a sloped region, a dead-band around nominal, a second sloped region,
and full absorption above the upper knee. See [`ieee1547_curve`](@ref) for the
curve used throughout the documentation.
"""
struct DroopCurve
    Vbp::Vector{Float64}
    qshape::Vector{Float64}
    function DroopCurve(Vbp, qshape)
        length(Vbp) == 6 || throw(ArgumentError("expected 6 breakpoint voltages, got $(length(Vbp))"))
        length(qshape) == 6 || throw(ArgumentError("expected 6 set-points, got $(length(qshape))"))
        issorted(Vbp) || throw(ArgumentError("breakpoint voltages must be non-decreasing"))
        new(collect(Float64, Vbp), collect(Float64, qshape))
    end
end

"""
    ieee1547_curve()

The Volt-VAr characteristic used throughout the documentation: breakpoint voltages
`[0.88, 0.90, 0.97, 1.00, 1.02, 1.10]` p.u. with shape `[1, 1, 0, 0, -1, -1]`, giving a
dead-band from 0.97 to 1.00 p.u. and saturation at ``\\pm\\bar q`` outside the knees.
"""
ieee1547_curve() = DroopCurve([0.88, 0.90, 0.97, 1.00, 1.02, 1.10], [1, 1, 0, 0, -1, -1])

"""
    droop_q(curve, v, qbar)

Evaluate the droop directly: the reactive output an inverter with capability `qbar`
would produce at terminal voltage `v`. Used to verify that an optimised dispatch
really does lie on the curve.
"""
function droop_q(curve::DroopCurve, v::Real, qbar::Real)
    V, q = curve.Vbp, curve.qshape
    v <= V[1] && return qbar * q[1]
    v >= V[6] && return qbar * q[6]
    for b in 1:5
        if v <= V[b+1]
            V[b+1] == V[b] && return qbar * q[b+1]
            θ = (v - V[b]) / (V[b+1] - V[b])
            return qbar * (q[b] + θ * (q[b+1] - q[b]))
        end
    end
    return qbar * q[6]
end

"""
    add_droop!(model, method, curve, v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar)

Attach the Volt-VAr law to `model`, tying each inverter's reactive output `Qdg[d,h,m]`
to its terminal voltage `v[d,h,m]`. `method` selects the encoding:

- `:bigm`: one binary per segment activates that segment's voltage window and its
  affine law. Mixed-integer linear.
- `:lambda`: the operating point is a convex combination of the breakpoints, with
  SOS2 logic forcing the two active weights to be adjacent. Mixed-integer linear.
- `:heaviside`: segment masks built from unit steps, summed into one closed-form
  algebraic expression. Integer-free, but non-smooth, so it needs an NLP solver.

All three describe exactly the same curve and, given the same host model, return the
same dispatch. They differ in the solver technology they require and in how they
scale with the number of inverters and time steps.
"""
add_droop!(model, method::Symbol, args...) = add_droop!(model, Val(method), args...)

# ---------------------------------------------------------------------- Big-M ------
function add_droop!(model, ::Val{:bigm}, curve::DroopCurve,
                    v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar; Mbig = 1.1)
    Vbp = curve.Vbp
    # δ[i] = 1 selects segment i. W2 := δ₂·v and W4 := δ₄·v are exact linearisations
    # of a binary-times-bounded-continuous product. The tightest valid big-M here is
    # the largest attainable voltage; a loose M only weakens the LP relaxation.
    @variable(model, δ[1:5, DG_SET, HOUR_SET, QUARTER_SET], Bin)
    @variable(model, W2[DG_SET, HOUR_SET, QUARTER_SET])
    @variable(model, W4[DG_SET, HOUR_SET, QUARTER_SET])

    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        sum(δ[i, d, h, m] for i in 1:5) == 1)

    # flat segments 1, 3 and 5: the binary only has to switch on a voltage window
    for (i, lo, hi) in ((1, 1, 2), (3, 3, 4), (5, 5, 6))
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            v[d, h, m] >= Vbp[lo] - Mbig * (1 - δ[i, d, h, m]))
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            v[d, h, m] <= Vbp[hi] + Mbig * (1 - δ[i, d, h, m]))
    end

    # sloped segments 2 and 4: W = δ·v, whose bounds double as the voltage window
    for (i, W, lo, hi) in ((2, W2, 2, 3), (4, W4, 4, 5))
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            v[d, h, m] - W[d, h, m] >= -Mbig * (1 - δ[i, d, h, m]))
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            v[d, h, m] - W[d, h, m] <=  Mbig * (1 - δ[i, d, h, m]))
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            W[d, h, m] >= Vbp[lo] * δ[i, d, h, m])
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            W[d, h, m] <= Vbp[hi] * δ[i, d, h, m])
    end

    # the affine law of whichever segment is active; every coefficient is constant
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Qdg[d, h, m] ==
            δ[1, d, h, m] * qbar[d]
          + W2[d, h, m] * (-qbar[d] / (Vbp[3] - Vbp[2]))
          + δ[2, d, h, m] * (qbar[d] * Vbp[3] / (Vbp[3] - Vbp[2]))
          + W4[d, h, m] * (-qbar[d] / (Vbp[5] - Vbp[4]))
          + δ[4, d, h, m] * (qbar[d] * Vbp[4] / (Vbp[5] - Vbp[4]))
          - δ[5, d, h, m] * qbar[d])
    return model
end

# -------------------------------------------------------------- Lambda / SOS2 ------
function add_droop!(model, ::Val{:lambda}, curve::DroopCurve,
                    v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar)
    Vbp = curve.Vbp
    qpts = Dict(d => qbar[d] .* curve.qshape for d in DG_SET)
    # v and Qdg share one set of weights λ, so the operating point is pinned to the
    # curve. The binaries z impose SOS2: at most two λ nonzero, and adjacent.
    @variable(model, λ[1:6, DG_SET, HOUR_SET, QUARTER_SET] >= 0)
    @variable(model, z[1:5, DG_SET, HOUR_SET, QUARTER_SET], Bin)

    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        sum(λ[i, d, h, m] for i in 1:6) == 1)
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        sum(z[i, d, h, m] for i in 1:5) == 1)
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        λ[1, d, h, m] <= z[1, d, h, m])
    @constraint(model, [i in 2:5, d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        λ[i, d, h, m] <= z[i-1, d, h, m] + z[i, d, h, m])
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        λ[6, d, h, m] <= z[5, d, h, m])

    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        v[d, h, m] == sum(λ[i, d, h, m] * Vbp[i] for i in 1:6))
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Qdg[d, h, m] == sum(λ[i, d, h, m] * qpts[d][i] for i in 1:6))
    return model
end

# ------------------------------------------------------------------ Heaviside ------
function add_droop!(model, ::Val{:heaviside}, curve::DroopCurve,
                    v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar)
    Vbp = curve.Vbp
    # Windows W_b = H(v - Vbp_b) - H(v - Vbp_{b+1}) are 1 on their own segment and 0
    # elsewhere, so the masked sum collapses to the single active law. The sloped
    # terms are anchored at their zero crossings. No extra variables at all.
    Hstep(x) = op_ifelse(op_greater_than_or_equal_to(x, 0), 1.0, 0.0)
    α1 = Dict(d => -qbar[d] / (Vbp[3] - Vbp[2]) for d in DG_SET)
    α2 = Dict(d => -qbar[d] / (Vbp[5] - Vbp[4]) for d in DG_SET)

    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        Qdg[d, h, m] ==
            qbar[d] * (Hstep(v[d,h,m] - Vbp[1]) - Hstep(v[d,h,m] - Vbp[2]))
          + α1[d] * (v[d,h,m] - Vbp[3]) * (Hstep(v[d,h,m] - Vbp[2]) - Hstep(v[d,h,m] - Vbp[3]))
          + α2[d] * (v[d,h,m] - Vbp[4]) * (Hstep(v[d,h,m] - Vbp[4]) - Hstep(v[d,h,m] - Vbp[5]))
          - qbar[d] * (Hstep(v[d,h,m] - Vbp[5]) - Hstep(v[d,h,m] - Vbp[6])))
    return model
end
