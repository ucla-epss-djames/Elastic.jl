module Elastic

using PhysicalConstants.CODATA2018: k_B
using BlockingMethod

include("atoms.jl")
include("static.jl")
include("voigt_royce.jl")
include("extrapolation.jl")
include("hashin_shtrikman.jl")
include("parrinello_rahman1981.jl")

export calculate_elastic_moduli

# right now designed to calculate the relevant strains for a cubic structure
function calculate_elastic_moduli(
    lattices::Matrix{Float64}, T::Union{Int64, Float64}
    )::Tuple{Matrix{Float64}, Matrix{Float64}}

    ϵ, V = strain_fluctuation(lattices)
    S = compliances_from_fluctuations(ϵ, V, T)
    Sc = symmetric_compliances(S)
    C = cubic_elastic_constants(Sc)
    mod = elastic_moduli(Sc, C, V, T)

    return (C, mod)
end

function elastic_moduli(S::Matrix{Float64}, C::Matrix{Float64}, V::Vector{Float64}, T::Union{Int64, Float64})::Matrix{Float64}

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

    # VOIGT-ROYCE G1 - DEBUGGING
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

    # VOIGT-ROYCE G2
    ABC = royce_bound(C[1,1], C[2,1], C[3,1])

    royce_G2 = mur(ABC[1], ABC[2], ABC[3])
    fx = -5.0 * (4.0*ABC[1] - 4.0*ABC[2] + 3.0*ABC[3])^(-2)
    royce_Gx2 = S[1,3]^2*(4.0*fx)^2 + S[2,3]^2*(-4.0*fx)^2 + S[3,3]^2*(3.0*fx)^2
    royce_Gx2 = sqrt(royce_Gx2)

    vr_G_avg2 = (voigt_G + royce_G2) / 2.0
    vr_Gx_avg2 = voigt_Gx^2 / 4.0 + royce_Gx2^2 / 4.0
    vr_Gx_avg2 = sqrt(vr_Gx_avg2)

    # HASHIN-SHTRIKMAN G - DOESNT WORK
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

    mod[1,:] = [voigt_K voigt_Kx]
    mod[2,:] = [K Kx]
    mod[2,:] = [vr_G_avg vr_Gx_avg]
    mod[3,:] = [vr_G_avg2 vr_Gx_avg2]
    mod[4,:] = [hs_G_avg hs_Gx_avg]

    return mod
end

end # end of module
