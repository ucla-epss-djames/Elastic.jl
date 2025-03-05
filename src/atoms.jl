using Printf

export nearest_neighbor_table


"Calculates the distance between two atoms of a crystal lattice."
function distance(r1, r0, lattice; periodic=false)
    # r1, r0 are 3 component vectors (x, y, z)
    # r is the distance between each point
    r = r1 .- r0
    # lattice is the components of the phase (a, b, c)
    if(periodic==true)
        r_scaled = (r .- round.(r)) .* lattice
    else
        r_scaled = r .* lattice
    end
    # dividing by two for the atom size on each end
    d = sqrt(sum(x -> x^2, r_scaled))
    return d
end

"Finds the nearest neighbor across all the atoms."
function nearest_neighbor_table(lattice, atoms; coord=0)

    if(coord != 0)
        atoms = coordination_shell(atoms, coord)
    end
    n = length(atoms)
    l = Int64(n*(n-1) / 2)
    table = zeros(l, 3)

    index = 1
    for i in 1:n-1
        for j in i+1:n
            val = distance(atoms[i][1], atoms[j][1], lattice)
            table[index,1] = val
            table[index,2] = i
            table[index,3] = j
            index += 1
        end
    end


    return (table, atoms)
end

function coordination_shell(atoms, coord)

    n = length(atoms)

    # supercell
    coord_range = -coord:coord
    coord_tuples = [(xi, yi, zi) for xi in coord_range, yi in coord_range, zi in coord_range]

    new_atoms = Dict{Int64, Any}()

    index = 1
    for (x, y, z) in coord_tuples
        for (atom_index, atom_data) in atoms
            new_coords = atom_data[1] .+ (x, y, z)
            new_atoms[index] = (new_coords, atom_data[2])
            index += 1
        end
    end

    return new_atoms
end

function display_table(table, atoms)

    l = length(table[:,1])

    for i in 1:l
        @printf("%3.5f %3s - %3s\n", table[i,1], atoms[table[i,2]][2], atoms[table[i,3]][2])
    end

end
