module Elastic

using LinearAlgebra

export c1x, c44, Kv, muv, Kr, mur, Cs, Gu, Gl, royce_bound

# Calculating Elastic Constants
c1x(uu::Real, xx::Real, e::Real) = (uu - xx) / e
c44(zx::Real, e::Real) = -zx / (2e)

# Voigt bound
Kv(c11::Real, c12::Real) = 1/3*(c11 + 2c12)
muv(c11::Real, c12::Real, c44::Real) = 1/5*(c11 - c12 + 3c44)

# Royce bound
Kr(a::Real, b::Real) = (3(a + 2b))^-1
mur(a::Real, b::Real, c::Real) = 5 / (4a - 4b + 3c)

# s constants for Royce bound calculations
ar(s11::Real, s22::Real, s33::Real) = 1/3*(s11 + s22 + s33)
br(s12::Real, s23::Real, s31::Real) = 1/3*(s12 + s23 + s31)
cr(s44::Real, s55::Real, s66::Real) = 1/3*(s44 + s55 + s66)

# Hashin and Shtrikman bounds
Cs(c11::Real, c12::Real) = (c11 - c12) / 2
Gu(K::Real, cs::Real, c44::Real) = c44 + 2*(5/(cs - c44) + 18*(K + 2c44)/(5c44*(3K + 4c44)))
Gl(K::Real, cs::Real, c44::Real) = cs + 3*(5/(c44 - cs) + 12*(K + 2cs)/(5cs*(3K + 4cs)))

function royce_bound(c11::Real, c12::Real, c44::Real)

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

end # end of module
