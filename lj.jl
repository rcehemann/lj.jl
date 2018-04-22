# run script for lj.jl
#----------------------------------------------------------------------------#
# settings are defined in the beginning section, simulation is run in the end
# ############################################################################

include("./src/sim.jl")
include("./src/pot.jl")
include("./src/box.jl")

#-----------------------------------INPUTS-----------------------------------#
# pot parameters
ljSigma = 1./1.122
ljEpsilon = 1.
ljCut = 4.
pot = Pot(ljSigma, ljEpsilon, ljCut)

# cell parameters
lat = "hex"
xDim = 10
yDim = 10
box = Box(lat, [xDim, yDim], pot)

# sim parameters
dt = 1e-6
nsteps = 1000
sim = Sim(box, dt, nsteps)
#----------------------------------------------------------------------------#


#--------------------------------------RUN-----------------------------------#
runSim(sim)
