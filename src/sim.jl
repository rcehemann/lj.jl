# simulation type and methods
#----------------------------------------------------------------------------#
# timestep : simulation timestep
# time : current simulation time
# tmax : end simulation time
# box : simulation box
# ############################################################################

include("./box.jl")

type Sim
        timestep::Float64
        time::Float64
        tmax::Float64
        box::Box
end

# constructor taking box, timestep and number of steps
function Sim(box::Box, dt, nsteps)
        Sim(dt, 0.0, dt * nsteps, box)
end

# perform a Velocity Verlet step
function vvStep(sim::Sim)
        f_t = map(x -> x.f, sim.box.atoms)
        v_t = map(x -> x.v, sim.box.atoms)
        r_t = map(x -> x.r, sim.box.atoms)
        m   = map(x -> x.m, sim.box.atoms)

        # compute and update positions
        r_t .+= sim.timestep .* v_t + (sim.timestep)^2 .* f_t ./ (2 .* m)
        updatePosition(sim.box, r_t)

        # compute and update first half-step in velocity (LeapFrog)
        v_t  .+= f_t ./ (2 .* m)
        updateVelocity(sim.box, v_t)

        # compute and update forces
        f_t = updateForce(sim.box, totalForce(sim.box.pot, sim.box.atoms))

        # compute and update second LeapFrog step in velocity
        v_t  .+= f_t ./ (2 .* m)
        updateVelocity(sim.box, v_t)

        # update time and energy
        updateEnergy(sim.box)
        sim.time += sim.timestep
end

# run the simulation, performing VV steps and printing total energies
function runSim(sim::Sim)
        while sim.time <= sim.tmax
            vvStep(sim)
            @printf "time: %.3f \t energy: %.12f\n" sim.time sim.box.energy
        end
        println("All done!")
end
