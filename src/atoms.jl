export nearest_neighbor_table

"Calculates the distance between two atoms of a crystal lattice."
function distance(r1, r0, lattice)
    # r1, r0 are 3 component vectors (x, y, z)
    # r is the distance between each point
    r = r1 .- r0
    # lattice is the components of the phase (a, b, c)
    r_scaled = (r .- round.(r)) .* lattice
    # dividing by two for the atom size on each end
    d = sqrt(sum(x -> x^2, r_scaled))
    return d
end

"Finds the nearest neighbor across all the atoms."
function nearest_neighbor_table(lattice, atoms)
    n = length(atoms)
    l = Int64(n*(n-1) / 2)
    d = zeros(l, 3)

    index = 1
    for i in 1:n-1
        for j in i+1:n
            val = distance(atoms[i], atoms[j], lattice)
            d[index,1] = val
            d[index,2] = i
            d[index,3] = j
            index += 1
        end
    end

    return d
end
