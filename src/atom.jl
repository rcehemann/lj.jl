# Atom type
#----------------------------------------------------------------------------#
# r : position vector
# v : velocity vector
# f : force vector
# m : mass
# t : atom species (presently unused)
# ############################################################################

type AtomLJ
    r::Array{Float64}
    v::Array{Float64}
    f::Array{Float64}
    m::Float64 # mass
    t::Int64  # type
end

# default constructor
function AtomLJ()
    AtomLJ(
            [0., 0.],
            [0., 0.],
            [0., 0.],
            12., 1
    )
end
