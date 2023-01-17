using Plots
using Trixi
using OrdinaryDiffEq

module TrixiExtensionRefine

using Trixi

struct IndicatorAlwaysRefine{Cache<:NamedTuple} <: Trixi.AbstractIndicator
  cache::Cache
end

function IndicatorAlwaysRefine(semi)
  basis = semi.solver.basis
  alpha = Vector{real(basis)}()
  cache = (; semi.mesh, alpha)

  return IndicatorAlwaysRefine{typeof(cache)}(cache)
end

function (indicator::IndicatorAlwaysRefine)(u::AbstractArray{<:Any,4},
                                            mesh, equations, dg, cache;
                                            t, kwargs...)
  alpha = indicator.cache.alpha
  resize!(alpha, nelements(dg, cache))

  alpha .= 0.0

  for element in 1:length(alpha)
    # Calculate domain center of element.
    center = 0.5*(cache.elements.node_coordinates[:,   1,   1, element]
                + cache.elements.node_coordinates[:, end, end, element])

    alpha[element] = center[1] < 0.0 ? 1.0 : 0.0
  end

  return alpha
end

end # module TrixiExtensionRefine

import .TrixiExtensionRefine

###############################################################################
# semidiscretization of the linear advection equation

polydeg = 3
inilevel = 3
maxlevel = 5
trees_per_dimension = 1 .* (1, 1)

tspan = (0.0, 2.0)

advection_velocity = (0.5, 0.5)

equations = LinearScalarAdvectionEquation2D(advection_velocity)

function my_initial_condition(x, t, equation::LinearScalarAdvectionEquation2D)

  # println("x = ", x)

  scalar = exp(-(x[1]^2 + x[2]^2)/0.01)
  # scalar = 1.0

  return SVector(scalar)
end

initial_condition = my_initial_condition

solver = DGSEM(polydeg=polydeg, surface_flux=flux_lax_friedrichs)

# Deformed rectangle that looks like a waving flag,
f1(s) = 0.5 .* SVector(-1.0, s)
f2(s) = 0.5 .* SVector( 1.0, s)
f3(s) = 0.5 .* SVector(s, -1.0 + 0.1 * sin( pi * s))
f4(s) = 0.5 .* SVector(s,  1.0 + 0.1 * sin( pi * s))
faces = (f1, f2, f3, f4)

# This creates a mapping that transforms [-1, 1]^2 to the domain with the faces defined above.
# It generally doesn't work for meshes loaded from mesh files because these can be meshes
# of arbitrary domains, but the mesh below is specifically built on the domain [-1, 1]^2.
Trixi.validate_faces(faces)
mapping = Trixi.transfinite_mapping(faces)

mesh = T8codeMesh(trees_per_dimension,polydeg=polydeg, initial_refinement_level=inilevel,
  mapping=mapping, periodicity=true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

amr_indicator = TrixiExtensionRefine.IndicatorAlwaysRefine(semi)

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

alive_callback = AliveCallback(analysis_interval=1)

# save_solution = SaveSolutionCallback(interval=100,
#                                      save_initial_solution=true,
#                                      save_final_solution=true,
#                                      solution_variables=cons2prim)

# amr_controller = ControllerThreeLevel(semi, amr_indicator,
#                                       base_level=inilevel,
#                                       med_level=99, med_threshold=99999,
#                                       max_level=inilevel+1, max_threshold=0.5)

amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable=first),
                                      base_level=inilevel,
                                      med_level=99, med_threshold=99999,
                                      max_level=maxlevel, max_threshold=0.1)

amr_callback = AMRCallback(semi, amr_controller,
                           interval=5,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

stepsize_callback = StepsizeCallback(cfl=0.8)

visualization_callback = VisualizationCallback(interval=10, clims=(0,1.1), show_mesh=true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        # save_solution,
                        amr_callback,
                        visualization_callback,
                        stepsize_callback);

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);

summary_callback() # print the timer summary
