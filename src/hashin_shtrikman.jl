export Cs, Gu, Gl

# *** Hashin and Shtrikman Bounds ***

# elastic constant calculation for HS bounds
function Cs(c11::Float64, c12::Float64)::Float64
    return (c11 - c12) / 2.0
end

function Gu(K::Float64, cs::Float64, c44::Float64)::Float64
    a = 5.0 / (cs - c44)
    b = (18.0 / 5.0*c44) * (K + 2.0*c44) / (3.0*K + 4.0*c44)
    return c44 + 2.0*(a + b)
end

function Gl(K::Float64, cs::Float64, c44::Float64)::Float64
    a = 5.0 / (c44 - cs)
    b = 12.0 * (K + 2.0*cs) / (5.0*cs * (3.0*K + 4.0*cs))
    return cs + 3.0*(a + b)
end
