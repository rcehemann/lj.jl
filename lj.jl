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
lat = "sq"
xDim = 10
yDim = 10
box = Box(lat, [xDim, yDim], pot)

# sim parameters
dt = 1e-4
nsteps = 100
sim = Sim(box, dt, nsteps)
#----------------------------------------------------------------------------#


#--------------------------------------RUN-----------------------------------#
runSim(sim)
pos = map(x->x.r, sim.box.atoms)
xpos = [pos[i][1] for i=1:length(pos)]
ypos = [pos[j][2] for j=1:length(pos)]
