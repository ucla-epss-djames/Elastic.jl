export G1HS, G2HS

# *** Hashin and Shtrikman Bounds ***

# elastic constant calculation for HS bounds

function G1(c11::Float64, c12::Float64)::Float64
    return 1.0/2.0 * (c11 - c12)
end
# G2 = C44

function βx(K::Float64, Gx::Float64)::Float64
    num = 3.0 * (K + 2.0*Gx)
    dem = 5.0*Gx * (3.0*K + 4.0*Gx)
    return - num / dem
end

function G1HS(K::Float64, c11::Float64, c12::Float64, c44::Float64)::Float64
    g1 = G1(c11, c12)
    g2 = c44
    β1 = βx(K, g1)

    a = 5.0 / (g2 - g1)
    b = 4.0 * β1
    return g1 + 3.0*(a - b)^(-1)

end

function G2HS(K::Float64, c11::Float64, c12::Float64, c44::Float64)::Float64
    g1 = G1(c11, c12)
    g2 = c44
    β2 = βx(K, g2)

    a = 5.0 / (g1 - g2)
    b = 6.0 * β2
    return g2 + 2.0 * (a - b)^(-1)
end
