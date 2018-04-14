# simulation type and methods
#----------------------------------------------------------------------------#
# timestep : simulation timestep
# time : current simulation time
# tmax : end simulation time
# box : simulation box
# ############################################################################

module Simulation

    using Box

    type Simulation
        timestep::Float
        time::Float
        tmax::Float
        box::Box
    end

    # constructor taking box, timestep and number of steps
    function Simulation(c::Box, dt, nsteps)
        Simulation(dt, 0.0, dt * nsteps, c)
    end

    # perform a Velocity Verlet step
    function vvStep(sim::Simulation)
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
        updateForce(sim.box, totalForce(sim.box.pot, sim.box.atoms))
        f_t = map(x -> x.f, sim.box.atoms)

        # compute and update second LeapFrog step in velocity
        v_t  .+= f_t ./ (2 .* m)
        updateVelocity(sim.box, v_t)
    end

    # run the simulation, performing VV steps and printing total energies
    function run(sim::Simulation)
        while sim.time <= sim.tmax
            vvStep(sim)
            sim.time += sim.timestep
            @printf "time: %.6f energy: %.6f" sim.time sim.box.energy
        end
        println("All done!")
    end

end # module
