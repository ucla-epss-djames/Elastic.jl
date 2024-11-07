using StatsBase
using LinearAlgebra
using BlockingMethod
using PhysicalConstants.CODATA2018: k_B

include("voigt_royce.jl")
include("hashin_shtrikman.jl")

export calculate_elastic_moduli


# right now designed to calculate the relevant strains for a cubic structure
function calculate_elastic_moduli(
    lattices::Matrix{Float64}, T::Union{Int64, Float64}
    )::Tuple{Matrix{Float64}, Matrix{Float64}}

    ϵ, V = strain_fluctuation(lattices)
    S = compliances_from_fluctuations(ϵ, V, T)
    Sc = symmetric_compliances(S)
    C = cubic_elastic_constants(Sc)
    mod = elastic_moduli(Sc, C)

    return (C, mod)
end

# lattice file, each row: l11, l12, l13, l22, l23, l33
function strain_fluctuation(
    lattices::Matrix{Float64}
    )::Tuple{Matrix{Float64}, Vector{Float64}}

    l = size(lattices)[1]
    h = UpperTriangular(zeros(3, 3))
    h0 = UpperTriangular(zeros(3, 3))
    ϵ = zeros(l,6)
    V = zeros(l)

    # beneath eq 2.27 parrinello+rahman1981
    # h0 and how its defined in hernandez+2001
    h0[1,1] = estimate(lattices[:,1])[1]
    h0[2,2] = estimate(lattices[:,4])[1]
    h0[3,3] = estimate(lattices[:,6])[1]
    h0[2,3] = estimate(lattices[:,5])[1]
    h0[1,3] = estimate(lattices[:,3])[1]
    h0[1,2] = estimate(lattices[:,2])[1]

    for i in 1:l
        # populating h with "strained" vol
        h[1,1] = lattices[i,1]; h[2,2] = lattices[i,4]; h[3,3] = lattices[i,6]
        h[2,3] = lattices[i,5]; h[1,3] = lattices[i,3]
        h[1,2] = lattices[i,2]

        G = transpose(h) * h

        # Volume
        # beneath eq. 4 hernandez2001
        V[i] = det(h)

        # eq 17 parrinello+rahman1982
        h0_inv = inv(h0)
        h0_inv_T = transpose(h0_inv)
        Gh0_inv = G * h0_inv
        ϵ0 = 0.5*(h0_inv_T*Gh0_inv - I(3))
        ϵ[i,1] = ϵ0[1,1]
        ϵ[i,2] = ϵ0[2,2]
        ϵ[i,3] = ϵ0[3,3]
        ϵ[i,4] = ϵ0[2,3]
        ϵ[i,5] = ϵ0[1,3]
        ϵ[i,6] = ϵ0[1,2]

    end

    return (ϵ, V)
end

function compliances_from_fluctuations(
    ϵ::Matrix{Float64}, V::Vector{Float64}, T::Union{Int64, Float64}
    )::Matrix{Float64}

    ϵ_avg = zeros(9,2)
    l = size(ϵ)[1]

    Am = 1e30
    GPa = 1e-9
    kB = k_B.val

    # computing <ϵ,ϵ>
    x_avg = zeros(l)
    for i in 1:9
        if(i <= 3)
            x_avg = ϵ[:,i] .* ϵ[:,i]
        elseif(i == 4)
            x_avg = ϵ[:,1] .* ϵ[:,2]
        elseif(i == 5)
            x_avg = ϵ[:,1] .* ϵ[:,3]
        elseif(i == 6)
            x_avg = ϵ[:,2] .* ϵ[:,3]
        else
            x_avg = ϵ[:,i-3] .* ϵ[:,i-3]
        end

        ϵ_avg[i,:] .= estimate(x_avg)
    end

    V_avg = estimate(V)
    f = V_avg[1] / (kB * T) / Am / GPa
    S = ϵ_avg .* f

    return S
end

function symmetric_compliances(S::Matrix{Float64})::Matrix{Float64}
    Sc = zeros(3,3)

    # computing symmetic S
    for i in 1:3
        j = 3*(i - 1) + 1
        Sc[i,1] = (S[j,1] + S[j+1,1] + S[j+2,1]) / 3.0
        Sc[i,2] = sqrt(S[j,2]^2 + S[j+1,2]^2 + S[j+2,2]^2) / 3.0
        Sc[i,3] = std([S[j,2], S[j+1,2], S[j+2,2]])
    end

    return Sc
end

# computes elastic constants using compliances
# CALCULATE THE SHEAR MOD
# returns a 4x4 matrix
# row: C11, C12, C44, K, G
# col value, std
function cubic_elastic_constants(S::Matrix{Float64})::Matrix{Float64}

    # computing elastic constants
    # col 1: value
    # col 2: std
    a = S[1,1] - S[2,1]
    b = S[1,1] + 2*S[2,1]

    C = zeros(3,3)
    # --- C11 ---
    C[1,1] = (S[1,1] + S[2,1]) / (a*b)
    C[1,2] = S[1,2]^2 * (-2/(3*a^2) - 1/(3*b^2))^2
    C[1,2] += S[2,2]^2 * (2/(3*a^2) - 2/(3*b^2))^2
    C[1,2] = sqrt(C[1,2])
    C[1,3] += S[2,3]^2 * (2/(3*a^2) - 2/(3*b^2))^2
    C[1,3] = sqrt(C[1,3])

    # --- C12 ---
    C[2,1] = -S[2,1] / (a*b)
    C[2,2] = S[1,2]^2 * (1/(3*a^2) - 1/(3*b^2))^2
    C[2,2] += S[2,2]^2 * (-1/(3*a^2) - 2/(3*b^2))^2
    C[2,2] = sqrt(C[2,2])
    C[2,3] += S[2,3]^2 * (-1/(3*a^2) - 2/(3*b^2))^2
    C[2,3] = sqrt(C[2,3])

    # --- C44 ---
    # Nye Textbook
    C[3,1] = 1 / (4*S[3,1])
    C[3,2] = S[3,2] / (4*S[3,1]^2)
    C[3,3] = S[3,3] / (4*S[3,1]^2)


    return C

end

function elastic_moduli(S::Matrix{Float64}, C::Matrix{Float64})::Matrix{Float64}

    mod = zeros(3,2)
    # VOIGT-ROYCE BOUND

    # --- K ---
    voigt_K = Kv(C[1,1], C[2,1])
    voigt_Kx = C[1,3]^2 * (1.0/3.0)^2 + C[2,3]^2 * (2.0/3.0)^2
    voigt_Kx = sqrt(voigt_Kx)

    royce_K = Kr(S[1,1], S[2,1])
    fx = 3.0*(S[1,1] + 2.0*S[2,1])
    royce_Kx = S[1,3]^2*(-3.0*fx)^2 + S[2,3]^2*(-6.0*fx)^2
    royce_Kx = sqrt(royce_Kx)

    vr_K_avg = (voigt_K + royce_K) / 2.0
    vr_Kx_avg = voigt_Kx^2 / 4.0 + royce_Kx^2 / 4.0
    vr_Kx_avg = sqrt(vr_Kx_avg)

    # --- G ---
    voigt_G = muv(C[1,1], C[2,1], C[3,1])
    voigt_Gx = C[1,3]^2*(0.2)^2 + C[2,3]^2*(-0.2)^2 + C[3,3]*(0.6)^2
    voigt_Gx = sqrt(voigt_Gx)

    royce_G = mur(S[1,1], S[2,1], S[3,1])
    fx = -5.0 * (4.0*S[1,1] - 4.0*S[2,1] + 3.0*S[3,1])^(-2)
    royce_Gx = S[1,3]^2*(4.0*fx)^2 + S[2,3]^2*(-4.0*fx)^2 + S[3,3]^2*(3.0*fx)^2
    royce_Gx = sqrt(royce_Gx)

    vr_G_avg = (voigt_G + royce_G) / 2.0
    vr_Gx_avg = voigt_Gx^2 / 4.0 + royce_Gx^2 / 4.0
    vr_Gx_avg = sqrt(vr_Gx_avg)

    # HASHIN-SHTRIKMAN G
    hs_cs = Cs(C[1,1], C[2,1])

    hs_Gu = Gu(voigt_K, hs_cs, C[3,1])
    hs_Gux = 0
    hs_Gux = sqrt(hs_Gux)

    hs_Gl = Gl(voigt_K, hs_cs, C[3,1])
    hs_Glx = 0
    hs_Glx = sqrt(hs_Glx)

    hs_G_avg = (hs_Gu + hs_Gl) / 2.0
    hs_Gx_avg = hs_Gux^2 / 4.0 + hs_Glx^2 / 4.0
    hs_G_avg = sqrt(hs_G_avg)

    mod[1,:] = [vr_K_avg vr_Kx_avg]
    mod[2,:] = [vr_G_avg vr_Gx_avg]
    mod[3,:] = [hs_G_avg hs_Gx_avg]

    return mod
end