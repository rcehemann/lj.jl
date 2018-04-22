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

type Box
    natoms::Int64
    atoms::Array{AtomLJ}
    lattice::Array{Float64}
    dimensions::Array{Int64}
    pot::Pot
    energy::Float64
end

# empty constructor
function Box()
    Box(0, [], [], [], Pot())
end

# constructor with lattice and dimensions
function Box(lat, dim::Array{Int64}, pot::Pot)
    c = Box()
    setLattice(c, lat)
    setDimensions(c, dim)
    setPot(pot)
    createAtoms(c)
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
    end
end

# setter for dimensions in lattice units
function setDimensions(box, dimensions)
    # default to square if no lattice exists
    if length(box.lattice) == 0
        setLattice(box, "sq")
        return setDimensions(box, dimensions)
    elseif typeof(dimensions) == Array
        box.dimensions = [dimensions[i]*box.lattice[i] for i=1:2]
    end
end

# default mass and type setter
function setMassesAndTypes(box)
    [atom.t = 1 for atom in box.atoms]
    [atom.m = 1 for atom in box.atoms]
end

# initialize atoms on the lattice of the box
function createAtoms(box)
    if length(box.dimensions) == 0
        set_dimensions(box, [2,2])
        return create_atoms(box)
    else
        box.natoms = prod(box.dimensions)
        [box.atoms[i+j].r = ([
                box.lattice[1] * i
                box.lattice[2] * j
        ]) for i=1:box.dimensions[1] for j=1:box.dimensions[2]]
    end
    setMassesAndTypes(box)
end

# update atom positions
function updatePosition(box, r::Array{Float64})
    [a[i].r = r[i] for i=1:size(box.atoms)]
end

# update atom velocities
function updateVelocity(box, v::Array{Float64})
    [a[i].v = v[i] for i=1:size(box.atoms)]
end

# update atom forces
function updateForce(box, f::Array{Float64})
    [a[i].f = f[i] for i=1:size(box.atoms)]
end

# compute and store the total energy
function updateEnergy(box)
    box.energy = totalEnergy(box.pot, box.atoms)
end
