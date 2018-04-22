# Potential type and methods
#----------------------------------------------------------------------------#
# params : LJ params sigma and epsilon
# cutoff : outer radial cutoff for short-circuiting calculations
# ecut   : pairwise potential is shifted down s.t. U(cutoff) = 0
# ############################################################################

include("./atom.jl")

type Pot
    params::Array{Float64}
    cutoff::Float64 # outer cutoff
    ecut::Float64   # energy at cutoff
end

# default parameter values
function Pot()
    Pot(1., 1., 4.) # defaults
end

# constructor with custom params -- compute ecut
function Pot(sigma, epsilon, cutoff)
    ecut = -1*(sigma/cutoff)^6
    ecut -= ecut^2
    ecut *= 4*epsilon

    Pot([sigma, epsilon], cutoff, ecut)
end

# return energy of pair for LJ
function pairwiseEnergy(p::Pot, ai::AtomLJ, aj::AtomLJ)
    rij = aj.r-ai.r
    rij = norm(rij)

    if rij > p.cutoff || rij < 1e-8 # no self-energy
        return 0.0
    end

    lj = -1*(p.params[1]/rij)^6
    lj -= lj^2
    lj *= 4*p.params[2]
    return lj - p.ecut # energy is zero at cutoff
end

# return force on atom i
function pairwiseForce(p::Pot, ai::AtomLJ, aj::AtomLJ)
    rij = aj.r - ai.r
    dir = rij/norm(rij)
    rij = norm(rij)

    if rij > p.cutoff || rij < 1e-8
        return [0.0, 0.0]
    end

    fij = (6/p.params[1]) * (p.params[1]/rij)^7 - (12/p.params[1]) * (p.params[1]/rij)^13
    fij *= -4 * p.params[2]
    return fij .* dir
end

# compute total energy of a set of atoms
function totalEnergy(p::Pot, atoms::Array{AtomLJ})
    0.5 * sum(map(
        a -> pairwiseEnergy(p, a[1], a[2]),
        [(ai, aj) for ai in atoms for aj in atoms]
    ))
end

# compute the net force on each atom
function totalForce(p::Pot, atoms::Array{AtomLJ})
    [0.5 .* sum(map(
        aj -> pairwiseForce(p, ai, aj), atoms))
    for ai in atoms]
end
