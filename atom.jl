# Atom type
#----------------------------------------------------------------------------#
# r : position vector
# v : velocity vector
# f : force vector
# m : mass
# t : atom species (presently unused)
# ############################################################################

type Atom
    r::Array{Float64}(2)
    v::Array{Float64}(2)
    f::Array{Float64}(2)
    m::Float64 # mass
    t::Int64  # type
end
