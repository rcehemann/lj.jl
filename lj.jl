# run script for lj.jl
#----------------------------------------------------------------------------#
# settings are defined in the beginning section, simulation is run in the end
# ############################################################################

include("sim.jl")
include("cell.jl")
include("pot.jl")

#-----------------------------------INPUTS-----------------------------------#
# pot parameters
ljSigma = 1
ljEpsilon = 1
ljCut = 4
pot = Potential(ljSigma, ljEpsilon, ljCut)

# cell parameters
lat = "sq"
xDim = 10
yDim = 10
cell = Cell(lat, [xDim, yDim], pot)

# sim parameters
dt = 0.001
nsteps = 1000
sim = Simulation(cell, dt, nsteps)
#----------------------------------------------------------------------------#


#--------------------------------------RUN-----------------------------------#
run(sim)
