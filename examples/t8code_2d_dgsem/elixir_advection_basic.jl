using Trixi
using OrdinaryDiffEq

polydeg = 3
inilevel = 5

equations = LinearScalarAdvectionEquation2D((0.5,0.5))

solver = DGSEM(polydeg=polydeg, surface_flux=flux_lax_friedrichs,volume_integral=VolumeIntegralWeakForm())

coordinates_min = (-0.5,-0.5)
coordinates_max = ( 0.5, 0.5)

mesh = T8codeMesh((1,1),polydeg=polydeg, initial_refinement_level=inilevel,
  coordinates_min=coordinates_min, coordinates_max=coordinates_max, periodicity=true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_convergence_test, solver)

ode = semidiscretize(semi, (0.0,1.0))

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
# save_solution = SaveSolutionCallback(interval=100,
#                                      solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=1.6)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
# callbacks = CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)
callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback)

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);

summary_callback()
