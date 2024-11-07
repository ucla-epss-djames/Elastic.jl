using LinearAlgebra

export Kv, muv, Kr, mur, royce_bound

# *** Voigt Bound ***

# voigt bulk modulus
function Kv(c11::Float64, c12::Float64)::Float64
    return 1.0/3.0 * (c11 + 2.0*c12)
end

# voigt shear modulus
function muv(c11::Float64, c12::Float64, c44::Float64)::Float64
    return 1.0/5.0 * (c11 - c12 + 3.0*c44)
end


# *** Royce Bound ***

# royce bulk modulus
function Kr(a::Float64, b::Float64)::Float64
    return (3.0 * (a + 2.0*b)) ^ -1
end

# royce shear modulus
function mur(a::Float64, b::Float64, c::Float64)::Float64
    return 5.0 / (4.0*a - 4.0*b + 3.0*c)
end

# compliance calculations for royce calculation
function ar(s11::Float64, s22::Float64, s33::Float64)::Float64
    return 1.0/3.0 * (s11 + s22 + s33)
end

function br(s12::Float64, s23::Float64, s31::Float64)::Float64
    return 1.0/3.0 * (s12 + s23 + s31)
end

function cr(s44::Float64, s55::Float64, s66::Float64)::Float64
    return 1.0/3.0 * (s44 + s55 + s66)
end

# calculates compliances for royce bulk and shear modulus
function royce_bound(c11::Float64, c12::Float64, c44::Float64)::Vector{Float64}

    C = zeros(6,6)


    # array of elastic constants
    C[1,1] = c11
    C[2,2] = c11
    C[3,3] = c11
    C[1,2] = c12
    C[1,3] = c12
    C[2,1] = c12
    C[2,3] = c12
    C[3,1] = c12
    C[3,2] = c12
    C[4,4] = c44
    C[5,5] = c44
    C[6,6] = c44

    S = zeros(6,6)
    try
        S = inv(C)
    catch e
        if isa(e, LinearAlgebra.SingularException)
        end
    end

    # Solving out Royce Bound
    a = ar(S[1,1], S[2,2], S[3,3])
    b = br(S[1,2], S[2,3], S[3,1])
    c = cr(S[4,4], S[5,5], S[6,6])

    return [a, b, c]

end
