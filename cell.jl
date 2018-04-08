# cell type and methods
#----------------------------------------------------------------------------#
# natoms : number of atoms
# atoms  : array of Atom types
# lattce : 2D lattice on which to create atoms
# dimensions : boundaries of cell in lattice units
# pot : Potential type for computing forces and energies
# energy : total energy of the cell
# ############################################################################

module Cell

    include("atom.jl")

    type Cell
        natoms::Int
        atoms::Array{Atom}
        lattice::Array{Float64}(2,2)
        dimensions::Array{Int64}(2)
        pot::Potential
        energy::Float64
    end

    # empty constructor
    function Cell()
        Cell(0, [], [], [], Potential())
    end

    # constructor with lattice and dimensions
    function Cell(lat, dim::Array{Int64}(2), pot::Potential)
        c = Cell()
        setLattice(c, lat)
        setDimensions(c, dim)
        setPotential(pot)
        createAtoms(c)
    end

    # setter for potentials
    function setPotential(cell, pot::Potential)
        cell.pot = pot
    end

    # setter for lattice
    function setLattice(cell, lattice)
        if typeof(lattice) == String
            if lattice == "sq"
                cell.lattice = [1.0 0.0; 0.0 1.0]
            elseif lattice == "hex"
                cell.lattice = [1.0 0.0; 0.5 sqrt(3)/2]
            end
        end
    end

    # setter for dimensions in lattice units
    function setDimensions(cell, dimensions)
        # default to square if no lattice exists
        if length(cell.lattice) == 0
            setLattice(cell, "sq")
            return setDimensions(cell, dimensions)
        elseif typeof(dimensions) == Array
            cell.dimensions = [dimensions[i]*cell.lattice[i] for i=1:2]
        end
    end

    # initialize atoms on the lattice of the cell
    function createAtoms(cell)
        if length(cell.dimensions) == 0
            set_dimensions(cell, [2,2])
            return create_atoms(cell)
        else
            cell.natoms = prod(cell.dimensions)
            [cell.atoms[i+j].r = ([
                    cell.lattice[1] * i
                    cell.lattice[2] * j
            ]) for i=1:cell.dimensions[1] for j=1:cell.dimensions[2]]
        end
    end

    # update atom positions
    function updatePosition(cell, r::Array{Float64})
        [a[i].r = r[i] for i=1:size(cell.atoms)]
    end

    # update atom velocities
    function updateVelocity(cell, v::Array{Float64})
        [a[i].v = v[i] for i=1:size(cell.atoms)]
    end

    # update atom forces
    function updateForce(cell, f::Array{Float64})
        [a[i].f = f[i] for i=1:size(cell.atoms)]
    end

    # compute and store the total energy
    function updateEnergy(cell)
        cell.energy = totalEnergy(cell.pot, cell.atoms)
    end

end # module
