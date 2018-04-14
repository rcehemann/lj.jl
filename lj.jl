# run script for lj.jl
#----------------------------------------------------------------------------#
# settings are defined in the beginning section, simulation is run in the end
# ############################################################################

include("./src/Box.jl")
include("./src/Pot.jl")
include("./src/Sim.jl")
using Box
using Pot
using Sim

#-----------------------------------INPUTS-----------------------------------#
# pot parameters
ljSigma = 1
ljEpsilon = 1
ljCut = 4
pot = Pot(ljSigma, ljEpsilon, ljCut)

# cell parameters
lat = "sq"
xDim = 10
yDim = 10
box = Box.Box(lat, [xDim, yDim], pot)

# sim parameters
dt = 0.001
nsteps = 1000
sim = Sim.Sim(box, dt, nsteps)
#----------------------------------------------------------------------------#


#--------------------------------------RUN-----------------------------------#
run(sim)
