using Dates
using SHA

const ESCAPE_VERSION::String = "1.0.0a"
const HEADER::String = """
  
     ▄████████    ▄████████  ▄████████    ▄████████    ▄███████▄    ▄████████ 
    ███    ███   ███    ███ ███    ███   ███    ███   ███    ███   ███    ███ 
    ███    █▀    ███    █▀  ███    █▀    ███    ███   ███    ███   ███    █▀  
   ▄███▄▄▄       ███        ███          ███    ███   ███    ███  ▄███▄▄▄     
  ▀▀███▀▀▀     ▀███████████ ███        ▀███████████ ▀█████████▀  ▀▀███▀▀▀     
    ███    █▄           ███ ███    █▄    ███    ███   ███          ███    █▄  
    ███    ███    ▄█    ███ ███    ███   ███    ███   ███          ███    ███ 
    ██████████  ▄████████▀  ████████▀    ███    █▀   ▄████▀        ██████████ 
"""

const DESCRIPTION::String = "A Julia code for computing Polanyi characteristic curves from adsorbent geometries."
const AUTHORS::String = "F. Stavarache, J. M. Vicent-Luna, S. Calero, H. Tanaka"


function compute_checksum(file::AbstractString)
    if isfile(file)
        return bytes2hex(sha1(read(file)))
    else
        return "File not found"
    end
end


function write_header(file::IO, input_file::AbstractString, framework_file::AbstractString, 
                      probe_file::AbstractString, ff_file::AbstractString)
    println(file, HEADER)
    println(file, "ESCAPE ", ESCAPE_VERSION, "running on Julia ", VERSION)
    println(file, DESCRIPTION)
    println(file, "Authors: ", AUTHORS)
    println(file, "")
    println(file, "Date and Time: ", now())
    println(file, "Host: ", gethostname())
    println(file, "Operating System: ", Sys.KERNEL)
    println(file, "CPU Architecture: ", Sys.ARCH)
    println(file, "")
    println(file, "Directory: ", pwd())
    println(file, "Input File Checksum: ", compute_checksum(input_file))
    println(file, "Framework Geometry Checksum: ", compute_checksum(framework_file))
    println(file, "Probe Geometry Checksum: ", compute_checksum(probe_file))
    println(file, "Force Field File Checksum: ", compute_checksum(ff_file))
end

function write_footer(file::IO)
    println(file, "")
    println(file, "Simulation is finished!")
    println(file, "Date: ", Dates.Date(Dates.now()))
    println(file, "Time: ", Dates.Time(Dates.now()))
end

function write_section(file::IO, title::String)
    message = rpad("==== $title ", 81, "=")
    println(file, message)
    println(file, "")
end

function write_subsection(file::IO, title::String)
    message = rpad("---- $title ", 81, "-")
    println(file, message)
    println(file, "")
end

function write_result(file::IO, title::String, result::Real)
    key = rpad(title, 40, " ")
    argument = lpad(result, 40, " ")
    println(file, key, " ", argument)
end

function write_xyz(file::IO, atoms::Vector{Atom})
    
    index = rpad("Index", 10, " ")
    species = rpad("Species", 19, " ")
    x = rpad("X [angstrom]", 16, " ")
    y = rpad("Y [angstrom]", 16, " ")
    z = rpad("Z [angstrom]", 16, " ")
    println(file, index, " ", species, " ", x, " ", y, " ", z)
    
    for (index, atom) in enumerate(atoms)
        index = rpad(index, 10, " ")
        species = rpad(atom.species, 19, " ")
        x = rpad(atom.x, 16, " ")
        y = rpad(atom.y, 16, " ")
        z = rpad(atom.z, 16, " ")
        println(file, index, " ", species, " ", x, " ", y, " ", z)
    end
    println(file, "")
end

function write_pbc(file::IO, A::Array{Float64})
    
    box = rpad("Box [i, j, k]", 30, " ")
    x = rpad("X [angstrom]", 16, " ")
    y = rpad("Y [angstrom]", 16, " ")
    z = rpad("Z [angstrom]", 16, " ")
    println(file, box, " ", x, " ", y, " ", z)
    
    for i in [-1, 0, 1], j in [-1, 0, 1], k in [-1, 0, 1]
        x, y, z = A * [i, j, k]
        
        i = rpad("$i", 3, " ")
        j = rpad("$j", 3, " ")
        k = rpad("$k", 24, " ")
        x = rpad("$(round(x, digits=8))", 16, " ")
        y = rpad("$(round(y, digits=8))", 16, " ")
        z = rpad("$(round(z, digits=8))", 16, " ")
        println(file, i, " ", j, " ", k, " ", x, " ", y, " ", z)
    end
    println(file, "")
end
