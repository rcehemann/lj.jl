# simulation type and methods
#----------------------------------------------------------------------------#
# timestep : simulation timestep
# time : current simulation time
# tmax : end simulation time
# cell : simulation cell
# ############################################################################

module Simulation

    include("cell.jl")

    type Simulation
        timestep::Float
        time::Float
        tmax::Float
        cell::Cell
    end

    # constructor taking cell, timestep and number of steps
    function Simulation(c::Cell, dt, nsteps)
        Simulation(dt, 0.0, dt * nsteps, c)
    end

    # perform a Velocity Verlet step
    function vvStep(sim::Simulation)
        f_t = map(x -> x.f, sim.cell.atoms)
        v_t = map(x -> x.v, sim.cell.atoms)
        r_t = map(x -> x.r, sim.cell.atoms)
        m   = map(x -> x.m, sim.cell.atoms)

        # compute and update positions
        r_t .+= sim.timestep .* v_t + (sim.timestep)^2 .* f_t ./ (2 .* m)
        updatePosition(sim.cell, r_t)

        # compute and update first half-step in velocity (LeapFrog)
        v_t  .+= f_t ./ (2 .* m)
        updateVelocity(sim.cell, v_t)

        # compute and update forces
        updateForce(sim.cell, totalForce(sim.cell.pot, sim.cell.atoms))
        f_t = map(x -> x.f, sim.cell.atoms)

        # compute and update second LeapFrog step in velocity
        v_t  .+= f_t ./ (2 .* m)
        updateVelocity(sim.cell, v_t)
    end

    # run the simulation, performing VV steps and printing total energies
    function run(sim::Simulation)
        while sim.time <= sim.tmax
            vvStep(sim)
            sim.time += sim.timestep
            @printf "time: %.6f energy: %.6f" sim.time sim.cell.energy
        end
        println("All done!")
    end

end # module
