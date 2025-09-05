# **** work from STAT MECH *****
# TODO: adapt to correlation and dos to module
# TODO: v0 fit for extrapolation schemes from LARS code
# TODO: entropy calculation from Lars Github
using Plots
using QuadGK
using Statistics
using Polynomials
using LinearAlgebra
using Interpolations
using DelimitedFiles

const RUN = "./runs/"
const RET = "./results/"

const npart = 125
const cycles = 3000
const k_B = 1.3806485e-26
const mass = 39.948 / (6.02214086e23)
const LJepsilon = 1.65e-24 / mass
const LJsigma = 0.34
const r_m = 2^(1/6) * LJsigma
const delt_t = 0.005

function main()

    # reduced_Temp = 1.1
    # reduced_density = 0.7
    # Temp = LJepsilon * reduced_Temp * mass / k_B
    box(rd) = (npart * LJsigma^3 / rd)^(1/3)

    p1 = (run=1, ρ=0.7, T=1.1, box=box(0.7), label="liquid 1")
    p2 = (run=2, ρ=0.85, T=1.1, box=box(0.85), label="liquid 2")
    p3 = (run=3, ρ=0.05, T=1.1, box=box(0.05), label="gas")
    p4 = (run=4, ρ=0.4, T=1.4, box=box(0.4), label="super critical fluid")
    p = (p1, p2, p3, p4)

    # quest_1(p)
    # quest_2a(p)
    # quest_2b(p)
    # quest_3(p)
    quest_4(p)

end

# velocity autocorrelation function
# Cᵥᵥ(t) = <v̂(0)⋅v̂(t)> / <v²>
function corr_func(A, t; n=1)

    N = size(A)[1]
    C = zeros(t)

    for dt in 0:t-1
        t0 = 1:n:(N - dt*n)
        t1 = t0 .+ dt*n

        for i in 1:n
            i == 1 ? ni = 0 : ni = 1
            t0 = t0 .+ ni
            t1 = t1 .+ ni

            C[dt+1] += sum(dot.(A[t0,:], A[t1,:]))
        end
    end

    return C
end

# MSD pg 249 Chandler
# R(t)²= <|r̂(t) - r̂(0)|²>
function MSD(A, t; n=1)

    N = size(A)[1]
    R = zeros(t)

    for dt in 0:t-1
        t0 = 1:n:N .- dt*n
        t1 = t0 .+ dt*n

        for i in 1:n
            i == 1 ? ni = 0 : ni = 1
            t0 = t0 .+ ni
            t1 = t1 .+ ni

            r0 = sqrt.(A[t0,1].^2 .+ A[t0,2].^2 .+ A[t0,3].^2)
            r1 = sqrt.(A[t1,1].^2 .+ A[t1,2].^2 .+ A[t1,3].^2)
            R[dt+1] += sum((r1 .- r0).^2)
        end
        R[dt+1] /= length(t0) * n
    end

    return R
end

# g(r) pg 199 Chandler
# g(r) = <ρ(0)ρ(r)> / ρ²
function pair_dist(A, t, l; n=1, bins=100)

    B = Bool.(zeros(n))
    dl = l / bins

    h = zeros(bins)

    for dt in 1:t
        if(dt == 1)
            a = A[1:n*dt,:]
        else
            a = A[n*(dt-1)+1:n*dt,:]
        end
        B = Bool.(zeros(n))
        for ni in 1:n-1
            B[1] = 1

            r1_r0 = (a[.!B,:] .- a[B,:]) ./ l
            dr = r1_r0 .- round.(r1_r0)
            r = sqrt.(dr[:,1].^2 .+ dr[:,2].^2 .+ dr[:,3].^2) .* l

            bin = ceil.(Int, r ./ dl)

            for b in bin
                h[b] += 2
            end

            a = a[.!B,:]
            B = B[.!B]
        end
    end

    g = zeros(bins)
    V = l^3
    ρ = n / V
    K = 4*π/3 * ρ
    for i in 1:bins
        g[i] = h[i] / (n * t)
        r0 = (i - 1) * dl
        r1 = r0 + dl
        h_id = K * (r1^3 - r0^3)
        g[i] /= h_id
    end

    return g
end

function quest_1(p)

    data = gather_data(p[1].run)
    t = collect(keys(data))
    perm = sortperm(t)
    t = t[perm]
    t_max = length(t)

    C = zeros(t_max, 4)

    plt = plot(title="Velocity Autocorrelation Function",
               xaxis=("t [ps]", font(15)),
               yaxis=("Cᵥᵥ(t)", font(15)),
               legendfontsize=13,
               dpi=300)
    for i in 1:4
        println("COMPUTING RUN: ", i)

        vx = [data[ti][:,4] for ti in t]
        vx = vcat(vx...)
        vy = [data[ti][:,5] for ti in t]
        vy = vcat(vy...)
        vz = [data[ti][:,6] for ti in t]
        vz = vcat(vz...)

        C[:,i] = corr_func([vx vy vz], t_max, n=npart)
        C[:,i] /= C[1,i]

        plt = plot!(t, C[:,i], label=p[i].label)

        data = gather_data(p[i].run)

    end

    savefig(RET * "p1.png")

end

function quest_2a(p)

    data = gather_data(p[1].run)
    t = collect(keys(data))
    perm = sortperm(t)
    t = t[perm]
    t_max = length(t)

    R = zeros(t_max, 4)

    plts = Vector{Any}(nothing, 4)

    for i in 1:4
        println("COMPUTING RUN: ", i)

        rx = [data[ti][:,1] for ti in t]
        rx = vcat(rx...)
        ry = [data[ti][:,2] for ti in t]
        ry = vcat(ry...)
        rz = [data[ti][:,3] for ti in t]
        rz = vcat(rz...)

        R[:,i] = MSD([rx ry rz], t_max, n=npart)

        if(i == 4)
            t_half = 1500:t_max
            pfit = fit(t[t_half], R[t_half,i], 1)
        else
            pfit = fit(t, R[:,i], 1)
        end

        D = round(pfit[1] / 6, digits=5)

        plts[i] = plot(t, R[:,i], label=p[i].label)
        plts[i] = plot!(t, pfit.(t), linestyle=:dot, color=:black,
                        label="D = " * string(D))

        data = gather_data(p[i].run)
    end

    plt = plot(plts..., layout=4, legend=:topleft,
               xaxis="t [ps]", yaxis="R²(t)", dpi=300)
    savefig(RET * "p2a.png")

end

function quest_2b(p)

    bins = 200

    data = gather_data(p[1].run)
    t = collect(keys(data))
    perm = sortperm(t)
    t = t[perm]
    t_max = length(t)

    g = zeros(bins, 4)

    plts = Vector{Any}(nothing, 4)

    for i in 1:4
        println("COMPUTING RUN: ", i)

        rx = [data[ti][:,1] for ti in t]
        rx = vcat(rx...)
        ry = [data[ti][:,2] for ti in t]
        ry = vcat(ry...)
        rz = [data[ti][:,3] for ti in t]
        rz = vcat(rz...)

        len = p[i].box / LJsigma
        g[:,i] = pair_dist([rx ry rz], t_max, len, n=npart, bins=bins)

        l = range(0, len, length=bins)
        plts[i] = plot(l[1:100], g[1:100,i], label=p[i].label)

        data = gather_data(p[i].run)
    end

    plt = plot(plts..., layout=4, legend=:topright,
               xaxis="r [σ]", yaxis="g(r)", dpi=300)
    savefig(RET * "p2b.png")
end

# force autocorrelation function
# C(t) = <f̂(0)⋅f̂(t)> / <f²>
function quest_3(p)

    data = gather_data(p[1].run)
    t = collect(keys(data))
    perm = sortperm(t)
    t = t[perm]
    t_max = length(t)

    C = zeros(t_max, 4)

    plt = plot(title="Force Autocorrelation Function",
               xaxis=("t [ps]", font(15)),
               yaxis=("C(t)", font(15)),
               legendfontsize=13,
               dpi=300)
    for i in 1:4
        println("COMPUTING RUN: ", i)

        ax = [data[ti][:,7] for ti in t]
        ax = vcat(ax...)
        ay = [data[ti][:,8] for ti in t]
        ay = vcat(ay...)
        az = [data[ti][:,9] for ti in t]
        az = vcat(az...)

        C[:,i] = corr_func([ax ay az] .* mass, t_max, n=npart)
        C[:,i] /= C[1,i]

        plt = plot!(t, C[:,i], label=p[i].label)

        data = gather_data(p[i].run)

    end

    savefig(RET * "p3.png")

end

function quest_4(p)

    data = gather_data(p[1].run)
    t = collect(keys(data))
    perm = sortperm(t)
    t = t[perm]
    t_max = length(t)

    C = zeros(t_max, 4)
    tlin = range(0, 15, length=t_max)

    for i in 1:4
        println("COMPUTING RUN: ", i)

        vx = [data[ti][:,4] for ti in t]
        vx = vcat(vx...)
        vy = [data[ti][:,5] for ti in t]
        vy = vcat(vy...)
        vz = [data[ti][:,6] for ti in t]
        vz = vcat(vz...)

        C[:,i] = corr_func([vx vy vz], t_max, n=npart)
        C[:,i] /= C[1,i]

        c_int = LinearInterpolation(tlin, C[:,i])
        println(quadgk(c_int, 0, 15)[1] / 3)

        data = gather_data(p[i].run)

    end

end

function gather_data(run)

    steps = readdlm(RUN * "steps." * string(run) * ".out")
    cells = readdlm(RUN * "cells." * string(run) * ".out")

    m = size(cells)[1]
    l = length(steps)

    d = Dict{Float64,Matrix{Float64}}()

    i = 0
    for t in steps

        i0 = 1 + 125*i
        i1 = i0 + 124
        d[t] = cells[i0:i1,1:9]
        i += 1

    end

    return d

end

main()
