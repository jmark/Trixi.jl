using Trixi
using OrdinaryDiffEq

polydeg = 3
inilevel = 0
trees_per_dimension = (1,1)

equations = LinearScalarAdvectionEquation2D((1.0,0.0))

solver = DGSEM(polydeg=polydeg, surface_flux=flux_lax_friedrichs,volume_integral=VolumeIntegralWeakForm())

# # Deformed rectangle that looks like a waving flag,
scale = 0.1

f1(s) = SVector(-1.0, s)
f2(s) = SVector( 1.0, s)
f3(s) = SVector(s, -1.0 + scale * sin( pi * s))
f4(s) = SVector(s,  1.0 + scale *  sin( pi * s))

# Simple scaling of the map.
# xscale = 1.0
# yscale = 1.0
# f1(s) = xscale * SVector(-1.0, s)
# f2(s) = xscale * SVector( 1.0, s)
# f3(s) = yscale * SVector(s, -1.0)
# f4(s) = yscale * SVector(s,  1.0)
# faces = (f1, f2, f3, f4)

faces = (f1, f2, f3, f4)
mapping = Trixi.transfinite_mapping(faces)

mesh = T8codeMesh(trees_per_dimension,polydeg=polydeg, initial_refinement_level=inilevel,
  mapping=mapping, periodicity=true)


function my_initial_condition(x, t, equation::LinearScalarAdvectionEquation2D)

  # println("x = ", x)

  # scalar = exp(-(x[1]^2 + x[2]^2)/0.2^2)
  # scalar = 1.0
  scalar = x[1]

  return SVector(scalar)
end

initial_condition = my_initial_condition

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

ode = semidiscretize(semi, (0.0,0.0))

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
# callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback)
callbacks = CallbackSet(stepsize_callback)

# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
#             dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep=false, callback=callbacks);

sol = solve(ode, Euler(),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);



# summary_callback()
