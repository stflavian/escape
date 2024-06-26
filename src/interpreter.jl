using Plots
using ProgressBars

include("constants.jl")

"""
Structure for storing atomic Lennard-Jones parameters and charge.
"""
struct AtomProperties
    epsilon::Real
    sigma::Real
    charge::Real
    mass::Real
end

"""
Structure for storing the framework's lattice parameters.
"""
struct FrameworkProperties
    a::Real
    b::Real
    c::Real
    alpha::Real
    beta::Real
    gamma::Real
end

"""
Structure for storing the Cartesian coordinates of the framework's atoms. 
"""
struct Atom
    species::String
    x::Real
    y::Real
    z::Real
end


"""
    lennard_jones_energy(sigma::Real, epsilon::Real, distance::Real)

Compute the 6-12 Lennard-Jones energy between two particles

# Arguments
- `sigma::Real`: the size parameter, usually the sum of the particle's radii.
- `epsilon::Real`: the depth of the potential well.
- `distance::Real`: the distance between the center of masses of the two particles.
"""
function lennard_jones_energy(sigma::Real, epsilon::Real, distance::Real)
    frac = sigma / distance    
    return 4 * epsilon * KB * (frac^12 - frac^6)
end


"""
    coloumb_energy(charge::Real, distance::Real)

Compute the electrostatic energy between two particles.

# Arguments
- `charge::Real`: the product of the charges of the two particles.
- `distance::Real`: the distance between the center of masses of the two particles.
"""
function coloumb_energy(charge::Real, distance::Real)
    return Q^2 * charge * KE * 1e10 / distance
end


"""
    read_input_file(path::String)

Read the input file and store the data in appropriate variables.

# Arguments
- `path::String`: the location of the input file.
"""
function read_input_file(path::String)

    atom_prop = Dict{SubString{String}, AtomProperties}()
    frame_prop = FrameworkProperties
    framework = Vector{Atom}()
    probe = SubString{String}
    
    input_file = readlines(path)

    for line in input_file
        
        # Check if line is empty
        if length(line) == 0
            continue
        end

        line = split(line)
        if line[1] == "ATOMPROP"

            species = line[2]
            epsilon = parse(Float64, line[3])
            sigma = parse(Float64, line[4])
            charge = parse(Float64, line[5])
            mass = parse(Float64, line[6])

            atom_prop[species] = AtomProperties(epsilon, sigma, charge, mass)
        
        elseif line[1] == "FRAMEPROP"
            
            a = parse(Float64, line[2])
            b = parse(Float64, line[3])
            c = parse(Float64, line[4])
            alpha = parse(Float64, line[5])
            beta = parse(Float64, line[6])
            gamma = parse(Float64, line[7])
            
            frame_prop = FrameworkProperties(a, b, c, alpha, beta, gamma)
        
        elseif line[1] == "FRAMEWORK"

            species = line[2]
            x = parse(Float64, line[3])
            y = parse(Float64, line[4])
            z = parse(Float64, line[5])

            frame = Atom(species, x, y, z)
            push!(framework, frame)

        elseif line[1] == "PROBE"
            
            probe = line[2]

        end
    end
    return atom_prop, frame_prop, framework, probe
end


"""
    compute_potential(atom_prop::Dict{SubString{String}, AtomProperties}, 
    frame_prop::FrameworkProperties, framework::Vector{Atom}, probe::SubString{String}, 
    size::Integer)

Compute the potential 

# Arguments
- `atom_prop::Dict{SubString{String}, AtomProperties}`: 
- `frame_prop::FrameworkProperties`:
- `framework::Vector{Atom}`:
- `probe::SubString{String}`:
- `size::Integer`: 

"""
function compute_potential(atom_prop::Dict{SubString{String}, AtomProperties}, 
    frame_prop::FrameworkProperties, framework::Vector{Atom}, probe::SubString{String}, 
    size::Integer)
    
    # Compute the transformation matrix for fractional to Cartesian coordinates
    a = frame_prop.a
    b = frame_prop.b
    c = frame_prop.c

    alpha = frame_prop.alpha * pi / 180
    beta = frame_prop.beta * pi / 180
    gamma = frame_prop.gamma * pi / 180

    alphastar = acos((cos(beta) * cos(gamma) - cos(alpha)) / sin(beta) / sin(gamma)) 
    A = [a  b * cos(gamma)  c * cos(beta);
        0  b * sin(gamma)  c * -1 * sin(beta) * cos(alphastar);
        0  0  c * sin(beta) * sin(alphastar)]
   
    # Compute the offset displacements needed for periodic boundary conditions
    pbc_offsets = []
    for i in [-1, 0, 1], j in [-1, 0, 1], k in [-1, 0, 1]
        offset = A * [i, j, k]
        push!(pbc_offsets, offset)
    end
    
    # Initialize arrays and assign parameters for probe
    s = range(start=0, stop=1, length=size)
    potential = zeros(size, size, size, 4)

    sig1 = atom_prop[probe].sigma
    eps1 = atom_prop[probe].epsilon
    q1 = atom_prop[probe].charge

    for (i, fa) in tqdm(enumerate(s)), (j, fb) in enumerate(s), (k, fc) in enumerate(s)
                
        coordinates = A * [fa, fb, fc]
        x = coordinates[1]
        y = coordinates[2]
        z = coordinates[3]
        
        potential[i, j, k, 1] = x
        potential[i, j, k, 2] = y
        potential[i, j, k, 3] = z

        for atom in framework
        
            sig2 = atom_prop[atom.species].sigma
            eps2 = atom_prop[atom.species].epsilon
            q2 = atom_prop[atom.species].charge
            
            # Lorentz-Berthelot mixing rules and charge product
            sig = (sig1 + sig2) / 2
            eps = sqrt(eps1 * eps2)
            q = q1 * q2 
            
            for offset in pbc_offsets

                # Compute the position of the atomic image using PBC
                pb_atom_x = atom.x + offset[1]
                pb_atom_y = atom.y + offset[2]
                pb_atom_z = atom.z + offset[3]

                r = sqrt((pb_atom_x - x)^2 + (pb_atom_y - y)^2 + (pb_atom_z - z)^2)
            
                if  0.5 * sig < r < 5 * sig
                    potential[i, j, k, 4] += lennard_jones_energy(sig, eps, r)
                    potential[i, j, k, 4] += coloumb_energy(q, r)
                elseif r < 0.5 * sig
                    potential[i, j, k, 4] = 0
                    @goto finish_potenial_calculation
                elseif r > 5 * sig
                    continue
                end

            end
        end

        if potential[i, j, k, 4] > 0
            potential[i, j, k, 4] = 0
        end
        
        @label finish_potenial_calculation

        # Conversion from J to kJ/mol
        potential[i, j, k, 4] *= NA * 1e-3
    end
   
    mkdir("Output")

    x = zeros(size^3)
    y = zeros(size^3)
    z = zeros(size^3)
    pot = zeros(size^3)
    for i in 1:1:size, j in 1:1:size, k in 1:1:size
        index = i + (j-1) * size + (k-1) * size^2
        x[index] = potential[i, j, k, 1]
        y[index] = potential[i, j, k, 2]
        z[index] = potential[i, j, k, 3]
        pot[index] = potential[i, j, k, 4]
    end

    p = scatter(x, y, z, marker_z=pot, aspect_ratio=:equal, markersize=2, camera=(0, -90))
    xlabel!(p, "X [\$\\AA\$]")
    ylabel!(p, "Y [\$\\AA\$]")
    zlabel!(p, "Z [\$\\AA\$]")
    
    savefig(p, "Output/potential_landscape.png")

    return potential
end


"""
asda
"""
function compute_characteristic(atom_prop::Dict{SubString{String}, AtomProperties}, 
    frame_prop::FrameworkProperties, framework::Vector{Atom}, potential::Array{Float64, 4}, 
    size::Integer) 

    npoints = 30

    # Compute unitcell mass in g
    unitcell_mass = 0.0
    for atom in framework
        unitcell_mass += atom_prop[atom.species].mass
    end

    unitcell_mass = unitcell_mass * MC * 1e3
   
    # Compute unit cell volume
    a = frame_prop.a
    b = frame_prop.b
    c = frame_prop.c

    alpha = frame_prop.alpha * pi / 180
    beta = frame_prop.beta * pi / 180
    gamma = frame_prop.gamma * pi / 180

    unitcell_volume = a * b * c * sqrt(sin(alpha)^2 + sin(beta)^2 + sin(gamma)^2 + 
            2 * cos(alpha) * cos(beta) * cos(gamma) - 2)
    
    # Compute the volume of a sample point in ml
    sample_volume = unitcell_volume * 1e-24 / size^3
    
    minimum_potential = minimum(potential)
    potential_range = range(start=minimum_potential, stop=-0.000001, length=npoints)
    
    output_file = open("Output/characteristic.dat", "w+")
    write(output_file, "# Potential [kJ/mol] \t Volume [ml/g] \n")

    measured_potential = zeros(npoints)
    measured_volumes = zeros(npoints)
    for (index, ads_potential) in enumerate(potential_range)
        counter = 0
        for i in 1:1:size, j in 1:1:size, k in 1:1:size
            if potential[i, j, k, 4] <= ads_potential
                counter += 1
            end
        end
        
        # Store the volume in ml/g
        volume = counter * sample_volume / unitcell_mass
        measured_volumes[index] = volume
        
        # Store the positive value of potential in kJ/mol
        measured_potential[index] = -ads_potential
        
        write(output_file, "$(-ads_potential) \t $volume \n") 
    end

    close(output_file)

    characteristic_plot = plot(measured_potential, measured_volumes)
    xlabel!(characteristic_plot, "Potential [kJ/mol]")
    ylabel!(characteristic_plot, "Volume [ml/g]")
    savefig(characteristic_plot, "Output/characteristic.png") 
end
