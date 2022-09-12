using Trixi
using OrdinaryDiffEq
using Printf

###############################################################################
# semidiscretization of the linear advection equation

inilevel = 2
maxlevel = 4
polydeg = 3
trees_per_dimension = 1 .* (1, 1)

tspan = (0.0, 2.0)
# tspan = (0.0, 2.3)
# tspan = (0.0, 0.15)
# tspan = (0.0, 0.0)

advection_velocity = (0.5, 0.5)
# advection_velocity = (0.5, 0.0)
equations = LinearScalarAdvectionEquation2D(advection_velocity)

function my_initial_condition(x, t, equation::LinearScalarAdvectionEquation2D)

  # println("x = ", x)

  # scalar = x[1] + x[2]
  # scalar = exp(-(x[1]^2)/0.1)
  scalar = exp(-(x[1]^2 + x[2]^2)/0.1)
  # scalar = 0.0
  # scalar = 1.0*exp(-((x[1]-0.5)^2 + (x[2]-0.5)^2)/0.01)
  # scalar += 1.0*exp(-((x[1]-0.25)^2 + (x[2])^2)/0.01)
  # scalar += 1.0*exp(-((x[1]+0.25)^2 + (x[2])^2)/0.01)

  # scalar += 2.0*exp(-((x[1]+0.25)^2 + (x[2]-0.25)^2)/0.01)
  # scalar += 3.0*exp(-((x[1]-0.25)^2 + (x[2]+0.25)^2)/0.01)
  # scalar += 4.0*exp(-((x[1]+0.25)^2 + (x[2]+0.25)^2)/0.01)

  # scalar = 1.0

  return SVector(scalar)
end

initial_condition = my_initial_condition

boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = Dict(
  :all => boundary_condition,
  # :x_neg => boundary_condition,
  # :x_pos => boundary_condition,
  # :y_neg => boundary_condition,
  # :y_pos => boundary_condition,
  # Symbol("---") => boundary_condition
)

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

# Deformed rectangle that looks like a waving flag,
f1(s) = SVector(-1.0 + 0.1 * sin( pi * s), s)
f2(s) = SVector( 1.0 + 0.1 * sin( pi * s), s)
f3(s) = SVector(s, -1.0 + 0.1 * sin( pi * s))
f4(s) = SVector(s,  1.0 + 0.1 * sin( pi * s))
faces = (f1, f2, f3, f4)

# This creates a mapping that transforms [-1, 1]^2 to the domain with the faces defined above.
# It generally doesn't work for meshes loaded from mesh files because these can be meshes
# of arbitrary domains, but the mesh below is specifically built on the domain [-1, 1]^2.
Trixi.validate_faces(faces)
mapping_flag = Trixi.transfinite_mapping(faces)

coordinates_min = (-1.0,  1.0)
coordinates_max = ( 1.0, -1.0)

mapping_flip = Trixi.coordinates2mapping(coordinates_min, coordinates_max)

mapping(x,y) = mapping_flag(mapping_flip(x,y)...)

# mapping_flag = nothing

# Unstructured mesh with 24 cells of the square domain [-1, 1]^n
# mesh_file = joinpath(@__DIR__, "square_unstructured_2.inp")
# isfile(mesh_file) || download("https://gist.githubusercontent.com/efaulhaber/63ff2ea224409e55ee8423b3a33e316a/raw/7db58af7446d1479753ae718930741c47a3b79b7/square_unstructured_2.inp", mesh_file)

mesh_file = joinpath(@__DIR__,"meshfiles/unstructured_quadrangle.msh")

mesh = T8codeMesh{2}(mesh_file, polydeg=polydeg,
                    mapping=mapping,
                    initial_refinement_level=inilevel)

# mesh = T8codeMesh(
#   trees_per_dimension,
#   polydeg=polydeg, initial_refinement_level=inilevel, mapping=mapping, periodicity=true)
#   # polydeg=polydeg, initial_refinement_level=inilevel, mapping=mapping, periodicity=false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_restart = SaveRestartCallback(interval=100,
                                   save_final_restart=true)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable=first),
                                      base_level=inilevel,
                                      med_level=99, med_threshold=99999,
                                      max_level=maxlevel, max_threshold=0.1)

amr_callback = AMRCallback(semi, amr_controller,
                           interval=5,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

stepsize_callback = StepsizeCallback(cfl=0.8)

function my_save_plot(plot_data, variable_names;
                   show_mesh=true, plot_arguments=Dict{Symbol,Any}(),
                   time=nothing, timestep=nothing)

  Plots.plot(plot_data,clim=(0,1.1),title="Advected Blob");
  Plots.plot!(getmesh(plot_data))

  # Determine filename and save plot
  mkpath("out")
  filename = joinpath("out", @sprintf("solution_%06d.png", timestep))
  Plots.savefig(filename)
end

visualization_callback = VisualizationCallback(plot_creator=my_save_plot,interval=10, clims=(0,1.1), show_mesh=true)

callbacks = CallbackSet(# summary_callback,
                        # analysis_callback,
                        alive_callback,
                        # save_restart,
                        # save_solution,
                        amr_callback,
                        visualization_callback,
                        stepsize_callback
);

###############################################################################
# Run the simulation.

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # Solve needs some value here but it will be overwritten by the stepsize_callback.
            save_everystep=false, callback=callbacks);

summary_callback()
