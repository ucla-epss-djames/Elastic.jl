module Elastic

using PhysicalConstants.CODATA2018: k_B
using BlockingMethod

include("atoms.jl")
include("static.jl")
include("voigt_reuss.jl")
include("extrapolation.jl")
include("hashin_shtrikman.jl")
include("parrinello_rahman1981.jl")

export calculate_elastic_moduli

# right now designed to calculate the relevant strains for a cubic structure
function calculate_elastic_moduli(
    lattices::AbstractMatrix{<:Float64}, T::Real
    )::Tuple{Matrix{Float64}, Matrix{Float64}}

    ϵ, V = strain_fluctuation(lattices)
    S = compliances_from_fluctuations(ϵ, V, T)
    Sc = symmetric_compliances(S)
    C = cubic_elastic_constants(Sc)
    mod = elastic_moduli(Sc, C, V, T)

    return (C, mod)
end

function elastic_moduli(S::AbstractMatrix{<:Float64},
                        C::AbstractMatrix{<:Float64},
                        V::AbstractVector{<:Float64}, T::Real
                        )::Matrix{Float64}

    mod = zeros(5,2)

    # --- K --- ELASTIC CONSTANTS
    voigt_K = Kv(C[1,1], C[2,1])
    voigt_Kx = C[1,3]^2 * (1.0/3.0)^2 + C[2,3]^2 * (2.0/3.0)^2
    voigt_Kx = sqrt(voigt_Kx)

    # --- K --- VOLUME FLUX
    Am = 1e30
    GPa = 1e-9
    kB = k_B.val
    V_avg = estimate(V)
    V2_avg = estimate(V .^ 2)
    ΔV = V2_avg[1] - V_avg[1]^2
    fx = kB * T * Am * GPa
    K = V_avg[1] / ΔV[1] * fx
    Kx = V_avg[2]^2*(fx*(V2_avg[1]+V_avg[1]^2)/(V2_avg[1] - V_avg[1]^2)^2)^2
    Kx += V2_avg[2]^2*(-V_avg[1]*fx / (V2_avg[1] - V_avg[1]^2)^2)^2
    Kx = sqrt(Kx)

    # VOIGT-REUSS G1 - DEBUGGING
    voigt_G = muv(C[1,1], C[2,1], C[3,1])
    voigt_Gx = C[1,3]^2*(0.2)^2 + C[2,3]^2*(-0.2)^2 + C[3,3]*(0.6)^2
    voigt_Gx = sqrt(voigt_Gx)

    reuss_G = mur(S[1,1], S[2,1], S[3,1])
    fx = -5.0 * (4.0*S[1,1] - 4.0*S[2,1] + 3.0*S[3,1])^(-2)
    reuss_Gx = S[1,3]^2*(4.0*fx)^2 + S[2,3]^2*(-4.0*fx)^2 + S[3,3]^2*(3.0*fx)^2
    reuss_Gx = sqrt(reuss_Gx)

    vr_G_avg = (voigt_G + reuss_G) / 2.0
    vr_Gx_avg = voigt_Gx^2 / 4.0 + reuss_Gx^2 / 4.0
    vr_Gx_avg = sqrt(vr_Gx_avg)

    # VOIGT-REUSS G2
    ABC = reuss_bound(C[1,1], C[2,1], C[3,1])

    reuss_G2 = mur(ABC[1], ABC[2], ABC[3])
    fx = -5.0 * (4.0*ABC[1] - 4.0*ABC[2] + 3.0*ABC[3])^(-2)
    reuss_Gx2 = S[1,3]^2*(4.0*fx)^2 + S[2,3]^2*(-4.0*fx)^2 + S[3,3]^2*(3.0*fx)^2
    reuss_Gx2 = sqrt(reuss_Gx2)

    vr_G_avg2 = (voigt_G + reuss_G2) / 2.0
    vr_Gx_avg2 = voigt_Gx^2 / 4.0 + reuss_Gx2^2 / 4.0
    vr_Gx_avg2 = sqrt(vr_Gx_avg2)

    # HASHIN-SHTRIKMAN G - DOESNT WORK
    hs_G1 = G1HS(voigt_K, C[1,1], C[2,1], C[3,1])
    hs_G1x = sqrt(0)

    hs_G2 = G2HS(voigt_K, C[1,1], C[2,1], C[3,1])
    hs_G2x = sqrt(0)

    hs_G_avg = (hs_G1 + hs_G2) / 2.0
    hs_Gx_avg = hs_G1x^2 / 4.0 + hs_G2x^2 / 4.0
    hs_G_avg = sqrt(hs_G_avg)

    mod[1,:] = [voigt_K voigt_Kx]
    mod[2,:] = [K Kx]
    mod[2,:] = [vr_G_avg vr_Gx_avg]
    mod[3,:] = [vr_G_avg2 vr_Gx_avg2]
    mod[4,:] = [hs_G_avg hs_Gx_avg]

    return mod
end

end # end of module
