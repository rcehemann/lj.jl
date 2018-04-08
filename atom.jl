# atom module            #
# author: rcehemann      #
# 4/1/2018               #
# ########################

type Atom
    r::Array{Float64}(2)
    v::Array{Float64}(2)
    f::Array{Float64}(2)
    m::Float64 # mass
    t::Int64   # type
end
