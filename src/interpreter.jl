# TODO: The functions are too big and are very messy. They should be split into smaller
# functions. The print statements should also be removed or contained, since now it makes
# reading very difficult.

using Plots

include("constants.jl")
include("force_fields.jl")
include("probe.jl")
include("framework.jl")
include("writer.jl")

"""
    compute_potential_landscape(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, probe::Probe, sizea::Integer, sizeb::Integer, sizec::Integer,
    cutoff::Real, save::SubString{String})

Compute the potential landscape inside the framework for the given probe molecule. 

The unit cell is divided in `sizea` x `sizeb` x `sizec` smaller boxes, at the center of
which the probe is inserted. The Lennard-Jones and Coloumb interactions between the 
probe and the atoms in the framework are then computed and added to the potential value
in that position. If the distance between the probe and one of the atoms in the 
framework is closer than `0.5 sigma`, the potential value is set to 1. If the potential
value in a box is positive, it is set to 0. If the molecule contains more than one atom
this procedure is repeated for different rotations of the molecule, and the results
are averaged. When taking into account the probe-framework interactions, periodic
boundary conditions are used. 

The results are stored in a 4-dimensional array, where the first 3 entries represent 
the x, y, ans z coordinates, and the last entry is the value of the potential

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
- `probe::Probe`: the vector containing the constituting atoms of the probe.
- `sizea::Integer`: the number of units in which the a lattice vector is divided. 
- `sizeb::Integer`: the number of units in which the b lattice vector is divided.
- `sizec::Integer`: the number of units in which the c lattice vector is divided.
- `cutoff::Real`: the potential cutoff used for the energy calculations.
- `rotations::Integer: the number of rotations used for the molecule.
- `output_file::IO: the file to which the results are written real-time.
- `save::SubString{String}`: "yes" if the potential in each point should be plotted and
saved.
"""
function compute_potential_landscape(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, probe::Probe, sizea::Int64, sizeb::Int64, sizec::Int64,
    cutoff::Float64, rotations::Int64, output_file::IO, save::SubString{String})
    
    write_section(output_file, "Energy landscape calculations")
   
    # Compute the offset displacements needed for periodic boundary conditions
    write_subsection(output_file, "Periodic boundary conditions")
    A = compute_conversion_matrix(framework)
    write_pbc(output_file, A) 
    
    write_subsection(output_file, "PBC Framework")
    pbc_framework = generate_pbc(framework)
    write_xyz(output_file, pbc_framework.atoms)
    
    # Initialize arrays and assign parameters for probe
    sx = range(start=0, step=1/sizea, length=sizea) .+ 1/(sizea * 2)
    sy = range(start=0, step=1/sizeb, length=sizeb) .+ 1/(sizeb * 2)
    sz = range(start=0, step=1/sizec, length=sizec) .+ 1/(sizec * 2)
    potential = zeros(sizea, sizeb, sizec)
     
    # Main loop
    total_stats = @timed for (index, _) in enumerate(1:rotations)
        
        write_subsection(output_file, "Run $index")
       
        angle = rand() * 2 * pi
        axis_i = 1 - 2 * rand()
        axis_j = 1 - 2 * rand()
        axis_k = 1 - 2 * rand()
        rotated_probe = rotate_probe(probe, angle, axis_i, axis_j, axis_k)
        
        write_subsection(output_file, "Rotated probe")
        write_xyz(output_file, rotated_probe.atoms)

        run_stats = @timed for pos in eachindex(IndexCartesian(), potential)
            
            # Convert fractional coordinates to Cartesian coordinates
            x, y, z = fractional_to_cartesian(A, sx[pos[1]], sy[pos[2]], sz[pos[3]])

            for pbc_atom in pbc_framework.atoms
                
                # Store framework atom properties
                sig2 = atom_properties[pbc_atom.species].sigma
                eps2 = atom_properties[pbc_atom.species].epsilon
                
                # Compute distance between probe center of mass and framework atom
                cmr = compute_distance(x - pbc_atom.x, y - pbc_atom.y, z - pbc_atom.z)
                
                # Check if the probe center of mass falls within the atomic radius
                # of a framework atom 
                if cmr <= 0.5 * sig2 
                    @inbounds potential[pos] = 1
                    @goto next_point
                end
                
                for probe_atom in rotated_probe.atoms
                    
                    # Store probe atom properties
                    sig1 = atom_properties[probe_atom.species].sigma
                    eps1 = atom_properties[probe_atom.species].epsilon
                    
                    # Lorentz-Berthelot mixing rules and charge product
                    sig = lorentz_berthelot_sigma(sig1, sig2)
                    eps = lorentz_berthelot_epsilon(eps1, eps2)

                    # Compute distance between probe atom and framework atom
                    r = compute_distance(probe_atom.x + x - pbc_atom.x,
                                         probe_atom.y + y - pbc_atom.y,
                                         probe_atom.z + z - pbc_atom.z)
                    
                    # Check if distance is longer than cutoff
                    if r < cutoff
                        @inbounds potential[pos] += lennard_jones_energy(sig, eps, r) - 
                                                    lennard_jones_energy(sig, eps, cutoff)
                    end
                end
            end
            @label next_point
        end
        
        write_subsection(output_file, "Run results")

        total_boxes = sizea * sizeb * sizec
        inaccessible_boxes = length(potential[potential .== 1])
        
        positive_boxes = length(potential[potential .> 0])
        total_pos_pot = sum(potential[potential .> 0]) * 1e-3 / index

        negative_boxes = length(potential[potential .<= 0])
        total_neg_pot = sum(potential[potential .<= 0]) * 1e-3 / index
        
        write_result(output_file, "Total boxes evaluated [count]", total_boxes)
        write_result(output_file, "Inaccessible box ratio [-]", inaccessible_boxes/total_boxes)
        write_result(output_file, "Negative potential box ratio [-]", negative_boxes/total_boxes)
        write_result(output_file, "Positive potential box ratio [-]", positive_boxes/total_boxes)
        write_result(output_file, "Average potential [kJ/mol]", (total_pos_pot + total_neg_pot) / total_boxes)
        write_result(output_file, "Average positive potential [kJ/mol]", total_pos_pot/positive_boxes)
        write_result(output_file, "Average negative potential [kJ/mol]", total_neg_pot/negative_boxes)
        println(output_file, " ")

        write_subsection(output_file, "Run performance statistics")
        write_result(output_file, "Time elapsed [s]", run_stats.time) 
        write_result(output_file, "GC time elapsed [s]", run_stats.gctime)
        write_result(output_file, "Memory allocated [bytes]", run_stats.bytes)
        println(output_file, " ")
    end
    
    write_section(output_file, "Total results")
    
    total_boxes = sizea * sizeb * sizec
    inaccessible_boxes = length(potential[potential .== 1])
        
    positive_boxes = length(potential[potential .> 0])
    total_pos_pot = sum(potential[potential .> 0]) * 1e-3 / rotations

    negative_boxes = length(potential[potential .<= 0])
    total_neg_pot = sum(potential[potential .<= 0]) * 1e-3 / rotations
    
    write_result(output_file, "Total boxes evaluated [count]", total_boxes)
    write_result(output_file, "Inaccessible box ratio [-]", inaccessible_boxes/total_boxes)
    write_result(output_file, "Negative potential box ratio [-]", negative_boxes/total_boxes)
    write_result(output_file, "Positive potential box ratio [-]", positive_boxes/total_boxes)
    write_result(output_file, "Average potential [kJ/mol]", (total_pos_pot + total_neg_pot) / total_boxes)
    write_result(output_file, "Average positive potential [kJ/mol]", total_pos_pot/positive_boxes)
    write_result(output_file, "Average negative potential [kJ/mol]", total_neg_pot/negative_boxes)
    println(output_file, " ")

    write_subsection(output_file, "Run performance statistics")
    write_result(output_file, "Time elapsed [s]", total_stats.time) 
    write_result(output_file, "GC time elapsed [s]", total_stats.gctime)
    write_result(output_file, "Memory allocated [bytes]", total_stats.bytes)
    println(output_file, " ")
    
    if save == "yes"
        
        mkpath("Output")
        
        num = sizea * sizeb
        x = zeros(num)
        y = zeros(num)
        pot = zeros(num)
        index = 1
        for pos in eachindex(IndexCartesian(), potential)
            
            if pos[3] != 1
                continue
            end

            x[index], y[index], _ = fractional_to_cartesian(A, sx[pos[1]], sy[pos[2]], sz[pos[3]])
            total_potential = potential[pos]
            
            if total_potential == 1
                pot[index] = 1
            elseif total_potential >= 0
                pot[index] = 0
            else
                pot[index] = potential[pos] * 10^-3 / rotations
            end
            index += 1
        end
        
        p = scatter(x[pot .< 0], y[pot .< 0], marker_z=pot[pot .< 0], markersize=0.8, markerstrokewidth=0, 
                    showaxis=false, right_margin=12Plots.mm, legend=false,
                    colorbar=true, c=:acton, clims=(-15, 0), grid=false, aspect_ratio=:equal)
        scatter!(x[pot .== 0], y[pot .== 0], color="#255E11", markersize=0.8, markerstrokewidth=0, z_order=:back)
        scatter!(x[pot .== 1], y[pot .== 1], color="#0D4C00", markersize=0.8, markerstrokewidth=0)

        savefig(p, "Output/potential_landscape.pdf")
    end
    return potential .* 1e-3 ./ rotations
end


"""
    compute_characteristic(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, potential::Array{Float64, 4}, sizea::Integer, sizeb::Integer,
    sizec::Integer, npoints::Integer, save::SubString{String})

Compute the characteristic curve from the potential landscape.

Create a range between the minimum potential value recorded and 0 containing `npoints`.
For each energy value in the range, count the number of boxes that have a lower 
potential energy value. Multiply the number of boxes by the volume of a box to obtain
the total volume enclosed by each equipotential line.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
- `potential::Array{Float64, 4}`: the 4-dimensional array containing the value of the
potential at different points in the framework.
- `sizea::Integer`: the number of units in which the a lattice vector is divided. 
- `sizeb::Integer`: the number of units in which the b lattice vector is divided.
- `sizec::Integer`: the number of units in which the c lattice vector is divided.
- `npoints::Integer`: the number of points used for the characteristic curve.
- `save::SubString{String}`: "yes" if the characteristic curve should be plotted and
saved.
"""
function compute_characteristic(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, potential::Array{Float64, 3}, sizea::Int64, sizeb::Int64,
    sizec::Int64, npoints::Int64, save::SubString{String}) 

    framework_mass = compute_framework_mass(atom_properties, framework)
    framework_volume = compute_framework_volume(framework)

    # Compute the volume of a sample point in ml
    sample_volume = framework_volume * 1e-24 / sizea / sizeb / sizec
    
    minimum_potential = minimum(potential)
    potential_range = range(start=minimum_potential, stop=0, length=npoints)
    
    mkpath("Output")
    output_file = open("Output/characteristic.dat", "w+")
    write(output_file, "# Potential [kJ/mol] \t Volume [ml/g] \n")

    measured_potential = zeros(npoints)
    measured_volumes = zeros(npoints)
    for (index, ads_potential) in enumerate(potential_range)
        
        counter = length(potential[potential .<= ads_potential])

        # Store the volume in ml/g
        volume = counter * sample_volume / framework_mass / 1000
        measured_volumes[index] = volume

        # Store the positive value of potential in kJ/mol
        measured_potential[index] = -ads_potential

        println(output_file, -ads_potential, "\t", volume) 
    end

    close(output_file)
    
    if save == "yes"
        characteristic_plot = scatter(measured_potential, measured_volumes, color="#0D4C00", 
                                      markersize=3, legend=false, grid=false)
        xlims!((0, maximum(measured_potential)))
        ylims!((0, maximum(measured_volumes) * 1.1))
        xlabel!("Potential [kJ/mol]")
        ylabel!("Volume [ml/g]")
        savefig(characteristic_plot, "Output/characteristic.pdf") 
    end
end
