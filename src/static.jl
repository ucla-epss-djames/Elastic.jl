export c1x, c44

# Calculating Elastic Constants with a strained and unstrained lattice
function c1x(uu::Float64, xx::Float64, e::Float64)::Float64
    return (uu - xx) / e
end

function c44(zx::Float64, e::Float64)::Float64
    return -zx / (2.0*e)
end
