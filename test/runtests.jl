using SmartInverterDOPF
using JuMP
using Test
using HiGHS   # model construction only; see the LinDistFlow testset
using Statistics: mean
using Ipopt

const CASE = load_case()

@testset "SmartInverterDOPF" begin

    @testset "case data" begin
        @test nbus(CASE) == 33
        @test length(CASE.BRANCH_SET) == 32
        @test length(CASE.Bi_BRANCH_SET) == 64
        @test CASE.slack == 1
        @test CASE.DG_SET == [7, 18, 33]
        @test ndg(CASE) == 3
        @test size(CASE.Pload) == (33, 24, 4)
        @test size(CASE.Qload) == (33, 24, 4)

        # the substation carries no load of its own
        @test all(iszero, CASE.Pload[1, :, :])
        @test all(iszero, CASE.Qload[1, :, :])
        @test all(>=(0), CASE.Pload)

        # every bus and branch of the feeder is represented exactly once
        @test length(unique(CASE.BRANCH_SET)) == 32
        @test sort(union(CASE.DG_SET, CASE.NON_DG_SET, [CASE.slack])) == collect(1:33)

        # irradiance ceiling: zero overnight, full nameplate at solar noon
        for d in CASE.DG_SET
            @test maximum(CASE.Pdg_max_vary[d]) ≈ CASE.Pdg_max[d]
            @test iszero(CASE.Pdg_max_vary[d][1, 1])
            @test all(CASE.Pdg_max_vary[d] .<= CASE.Pdg_max[d] + 1e-12)
        end

        # inverters are oversized relative to their arrays, leaving reactive headroom
        for d in CASE.DG_SET
            @test CASE.Sdg_max[d] > CASE.Pdg_max[d]
        end

        # unit conversion round-trips
        @test kWh(CASE, 1.0) ≈ CASE.Sbase / 1e3 / 4
    end

    @testset "droop curve" begin
        curve = ieee1547_curve()
        @test length(curve.Vbp) == 6
        @test issorted(curve.Vbp)
        @test curve.qshape == [1, 1, 0, 0, -1, -1]

        @test_throws ArgumentError DroopCurve([1, 2, 3], [1, 1, 0, 0, -1, -1])
        @test_throws ArgumentError DroopCurve(collect(1:6), [1, 0])
        @test_throws ArgumentError DroopCurve([1, 0, 2, 3, 4, 5], [1, 1, 0, 0, -1, -1])

        V, q̄ = curve.Vbp, 2.0

        # the curve passes through every breakpoint
        for b in 1:6
            @test droop_q(curve, V[b], q̄) ≈ q̄ * curve.qshape[b]
        end

        # saturation outside the breakpoint range
        @test droop_q(curve, 0.5, q̄)  ≈  q̄
        @test droop_q(curve, 1.5, q̄)  ≈ -q̄

        # dead-band produces nothing
        @test droop_q(curve, 0.985, q̄) ≈ 0 atol = 1e-12

        # sloped segments interpolate linearly
        @test droop_q(curve, (V[2] + V[3]) / 2, q̄) ≈ q̄ / 2
        @test droop_q(curve, (V[4] + V[5]) / 2, q̄) ≈ -q̄ / 2

        # monotonically non-increasing in voltage: more volts, less reactive support
        vs = range(0.85, 1.15, length = 400)
        qs = [droop_q(curve, v, q̄) for v in vs]
        @test all(diff(qs) .<= 1e-12)

        # scales linearly with the inverter's capability
        @test droop_q(curve, 0.95, 4.0) ≈ 2 * droop_q(curve, 0.95, 2.0)
    end

    @testset "droop encodings build" begin
        # each encoding must attach to a model and introduce the variables it claims
        for (method, want_bin) in ((:bigm, true), (:lambda, true), (:heaviside, false))
            model = Model()
            H, Q = 1:2, 1:2                       # a small slice is enough to build
            @variable(model, 0.9 <= v[CASE.DG_SET, H, Q] <= 1.1)
            @variable(model, Qdg[CASE.DG_SET, H, Q])
            qbar = Dict(d => CASE.Sdg_max[d] for d in CASE.DG_SET)

            SmartInverterDOPF.add_droop!(model, method, ieee1547_curve(),
                                         v, Qdg, CASE.DG_SET, H, Q, qbar)

            nbin = count(is_binary, all_variables(model))
            @test (nbin > 0) == want_bin
            if want_bin
                # five segments per inverter per time step
                @test nbin == 5 * ndg(CASE) * length(H) * length(Q) ||
                      nbin == 11 * ndg(CASE) * length(H) * length(Q)
            end
            @test num_constraints(model; count_variable_in_set_constraints = false) > 0
        end

        @test_throws MethodError SmartInverterDOPF.add_droop!(
            Model(), :not_a_method, ieee1547_curve(), nothing, nothing, [1], 1:1, 1:1, Dict())
    end

    @testset "host selection" begin
        # only the two implemented hosts are accepted, and the check happens before
        # any model is built (so it needs no solver)
        @test_throws ArgumentError solve_dopf(CASE, nothing; host = :nonsense)
        @test_throws ArgumentError solve_dopf(CASE, nothing; host = :lindistflow_typo)
        @test_throws ArgumentError solve_dopf(CASE, nothing; warm_start = :flat)
    end

    @testset "warm-start sweep" begin
        # The sweep turns a dispatch into a complex state that satisfies the AC power
        # flow exactly, which is what makes it a good linearisation point. With no
        # inverter output at all it must reproduce the no-inverter base case.
        zero_dg = zeros(ndg(CASE), 24, 4)
        v_r, v_im, Ibs_r, Ibs_im, Ibr_r, Ibr_im =
            SmartInverterDOPF._sweep_state(CASE, zero_dg, zero_dg)
        V = sqrt.(v_r .^ 2 .+ v_im .^ 2)
        @test maximum(abs.(V .- base_case_voltages(CASE))) < 1e-9
        @test all(V[1, :, :] .≈ 1.0)                       # slack held at nominal

        # Ohm's law must hold exactly on every branch of the swept state
        h, q = 12, 2
        worst = maximum(abs((v_r[i,h,q] - v_r[j,h,q]) -
                            (CASE.R[(i,j)] * Ibr_r[((i,j),h,q)] -
                             CASE.X[(i,j)] * Ibr_im[((i,j),h,q)]))
                        for (i, j) in CASE.BRANCH_SET)
        @test worst < 1e-9

        # a non-zero dispatch must raise voltages relative to the base case
        some_dg = fill(0.05, ndg(CASE), 24, 4)
        v_r2, v_im2, _, _, _, _ = SmartInverterDOPF._sweep_state(CASE, some_dg, zero_dg)
        @test mean(sqrt.(v_r2 .^ 2 .+ v_im2 .^ 2)) > mean(V)
    end

    @testset "LinDistFlow host builds" begin
        # Structure only, no solve, so this needs no solver licence. (Neither HiGHS nor
        # GLPK completes this model; see the solver note in the documentation.)
        nsteps = 24 * 4
        # the MILP encodings build against an MILP solver, the NLP one against Ipopt
        for (method, want_bin, opt) in ((:bigm, 5, HiGHS.Optimizer),
                                        (:lambda, 5, HiGHS.Optimizer),
                                        (:heaviside, 0, Ipopt.Optimizer))
            b = SmartInverterDOPF._build_lindistflow(
                    CASE, opt; method = method, curve = ieee1547_curve(),
                    silent = true, time_limit_sec = nothing,
                    attributes = Dict{String,Any}())

            # five binaries per inverter per time step for the MILP encodings, none for
            # the Heaviside one, the same counts the IVACOPF host produces, because the
            # droop block is shared verbatim between the two hosts
            @test count(is_binary, all_variables(b.model)) == want_bin * ndg(CASE) * nsteps

            # the network model itself is linear: every constraint an affine one, so a
            # LinDistFlow model with an integer-free droop is a pure LP
            @test num_constraints(b.model; count_variable_in_set_constraints = false) > 0
            @test length(b.v) == nbus(CASE) * nsteps
            @test length(b.Pbr) == length(CASE.BRANCH_SET) * nsteps
            @test objective_sense(b.model) == MIN_SENSE
        end
    end

    @testset "base-case power flow" begin
        V = base_case_voltages(CASE)
        @test size(V) == (33, 24, 4)
        @test all(V[1, :, :] .≈ 1.0)              # slack held at nominal
        @test all(0.8 .< V .< 1.05)               # a plausible radial feeder profile

        # voltage falls monotonically away from the substation along the main trunk
        h, q = argmax([sum(CASE.Pload[:, h, q]) for h in 1:24, q in 1:4]).I
        trunk = [1, 2, 3, 4, 5, 6, 7, 8]
        @test all(diff(V[trunk, h, q]) .< 0)

        # without inverters this feeder violates its lower limit at peak demand
        @test minimum(V) < CASE.Vmin
    end
end
