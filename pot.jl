# Potential type and methods
#----------------------------------------------------------------------------#
# params : LJ params sigma and epsilon
# cutoff : outer radial cutoff for short-circuiting calculations
# ecut   : pairwise potential is shifted down s.t. U(cutoff) = 0
# ############################################################################

module Pot

    include("cell.jl")

    type Potential
        params::Array{Float64}
        cutoff::Float64 # outer cutoff
        ecut::Float64   # energy at cutoff
    end

    # default parameter values
    function Potential()
        Potential(1, 1, 4) # defaults
    end

    # constructor with custom params -- compute ecut
    function Potential(sigma, epsilon, cutoff)
        ecut = -1*(sigma/cutoff)^6
        ecut -= ecut^2
        ecut *= 4*epsilon

        Potential([sigma, epsilon], cutoff, ecut)
    end

    # return energy of pair for LJ
    function pairwiseEnergy(p::Potential, ai::Atom, aj::Atom)
        rij = aj.r-ai.r
        rij = norm(rij)
        if rij > p.cutoff
            return 0.0
        end

        lj = -1*(params[1]/rij)^6
        lj -= lj^2
        lj *= 4*params[2]
        return lj - p.ecut # energy is zero at cutoff
    end

    # return force on atom i
    function pairwiseForce(p::Potential, ai::Atom, aj::Atom)
        rij = aj.r - ai.r
        dir = rij/norm(rij)
        rij = norm(rij)
        if rij > p.cutoff
            return [0.0, 0.0]
        end
        fij = 6 * (params[2]/rij)^7 - 12 * (params[1]/rij)^13
        fij *= -4 * epsilon
        return fij .* dir
    end

    # compute total energy of a set of atoms
    function totalEnergy(p::Potential, atoms::Array{Atom})
        0.5 * sum(map(
            (x, y) -> pairwiseEnergy(p, x, y), atoms, atoms
        ))
    end

    # compute the net force on each atom
    function totalForce(p::Potential, atoms::Array{Atom})
        0.5 .* [sum(map(
            (x, y) -> pairwiseForce(p, x, y), atom, atoms
        )) for atom in atoms]
    end

end # module
