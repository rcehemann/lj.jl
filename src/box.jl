# box type and methods
#----------------------------------------------------------------------------#
# natoms : number of atoms
# atoms  : array of Atom types
# lattce : 2D lattice on which to create atoms
# dimensions : boundaries of box in lattice units
# pot : Pot type for computing forces and energies
# energy : total energy of the box
# ############################################################################

include("./pot.jl")
include("./atom.jl")

type Box
    natoms::Int64
    atoms::Array{AtomLJ}
    lattice::Array{Float64,2}
    bounds::Array{Float64,2}
    dimensions::Array{Int64}
    pot::Pot
    energy::Float64
end

# empty constructor
function Box()
    Box(1, [AtomLJ()], [1. 0.; 0. 1.], [1. 0.; 0. 1.], [1, 1], Pot(), 0.)
end

# constructor with lattice and multiplicative dimensions
function Box(lat, dim::Array{Int64}, pot::Pot)
    box = Box()
    box.dimensions = dim
    setLattice(box, lat)
    setBounds(box, dim)
    setPotential(box, pot)
    createAtoms(box)

    return box # return box
end

# setter for Pots
function setPotential(box, pot::Pot)
    box.pot = pot
end

# setter for lattice
function setLattice(box, lattice)
    if typeof(lattice) == String
        if lattice == "sq"
            box.lattice = [1.0 0.0; 0.0 1.0]
        elseif lattice == "hex"
            box.lattice = [1.0 0.0; 0.5 sqrt(3)/2]
        end
    elseif typeof(lattice) == Array{Float64, 2}
        box.lattice = lattice
    end
end

# setter for dimensions in lattice units
function setBounds(box, dimensions)
    # default to square if no lattice exists
    if length(box.lattice) == 0
        setLattice(box, "sq")
        return setBounds(box, dimensions)
    else
        box.bounds = dimensions .* box.lattice
    end
end

# default mass and type setter
function setMassesAndTypes(box)
    [atom.t = 1  for atom in box.atoms]
    [atom.m = 10. for atom in box.atoms]
end

# initialize atoms on the lattice of the box
function createAtoms(box)
    if length(box.dimensions) == 0
        setBounds(box, [2,2])
        create_atoms(box)
    else
        box.natoms = prod(box.dimensions)
        box.atoms  = [AtomLJ() for i=1:box.natoms]
        updatePosition(box, [
                box.lattice[1,:] * i + box.lattice[2,:] * j
         for i=1:box.dimensions[1] for j=1:box.dimensions[2]])
    end
    setMassesAndTypes(box)
    box.energy = totalEnergy(box.pot, box.atoms)
end

# update atom positions
function updatePosition(box, r::Array{Array{Float64,1},1})
    [box.atoms[i].r = r[i] for i=1:box.natoms]
end

# update atom velocities
function updateVelocity(box, v::Array{Array{Float64,1},1})
    [box.atoms[i].v = v[i] for i=1:box.natoms]
end

# update atom forces
function updateForce(box, f::Array{Array{Float64,1},1})
    [box.atoms[i].f = f[i] for i=1:box.natoms]
end

# compute and store the total energy
function updateEnergy(box)
    box.energy = totalEnergy(box.pot, box.atoms)
end
