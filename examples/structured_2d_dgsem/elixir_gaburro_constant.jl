using Plots
using Trixi
using Printf
using OrdinaryDiffEq


equations = Gaburro2D(1.0, 2.78*10^5, 1000.0, 9.81)

function initial_condition_const(x, t, equations::Gaburro2D)
  # liquid domain
  rho = 1200.0
  v1 = 0.0
  v2 = 0.0
  alpha = 1.0
    
  return prim2cons(SVector(rho, v1, v2, alpha), equations)
end
  
initial_condition = initial_condition_const

function boundary_condition_wallG(u_inner, normal_direction::AbstractVector, direction, x, t,
                                  surface_flux_function, equations::Gaburro2D)
  
  if direction in (1, 2) # x direction
    u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])
  else # y direction
    u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4])
  end

  # println(surface_flux_function)
  surf_flux = surface_flux_function[1]

  # Calculate boundary flux
  if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
    flux = surf_flux(u_inner, u_boundary, normal_direction, equations)
  else # u_boundary is "left" of boundary, u_inner is "right" of boundary
    flux = surf_flux(u_boundary, u_inner, normal_direction, equations)
  end
  
  return flux
end
  

boundary_conditions = (x_neg=boundary_condition_wallG,
                       x_pos=boundary_condition_wallG,
                       y_neg=boundary_condition_wallG,
                       y_pos=boundary_condition_wallG)
  
# #volume_flux = (flux_central, flux_nonconservative_gaburro)
# solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)#(flux_lax_friedrichs, flux_nonconservative_gaburro),
#                  #volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

surface_flux = (flux_lax_friedrichs, flux_nonconservative_gaburro)
volume_flux = (flux_central, flux_nonconservative_gaburro)

solver = DGSEM(polydeg=1,
  surface_flux=surface_flux, 
  volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-0.5, 0.0) # minimum coordinates (min(x), min(y))
coordinates_max = ( 0.5, 1.0) # maximum coordinates (max(x), max(y))

cells_per_dimension = (16, 16)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity=false)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions,
                                    source_terms=source_terms_gravity)


tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

stepsize_callback = StepsizeCallback(cfl=0.1)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

function my_save_plot(plot_data, variable_names;
                   show_mesh=true, plot_arguments=Dict{Symbol,Any}(),
                   time=nothing, timestep=nothing)

  title = @sprintf("2D KHI | Trixi.jl | 2nd-order DG | Gaburro: t = %3.2f", time)

  sol = plot_data["rho"]

  Plots.plot(sol,
    clim=(0.0,1_500),
    colorbar_title="\ndensity",
    title=title,titlefontsize=10,
    dpi=300,
  )
  Plots.plot!(getmesh(plot_data),linewidth=0.5)

  mkpath("out")
  filename = joinpath("out", @sprintf("solution_%06d.png", timestep))
  Plots.savefig(filename)
end

visualization_callback = VisualizationCallback(
  plot_creator=my_save_plot,interval=1000, clims=(0,1.1), show_mesh=true)

callbacks = CallbackSet(stepsize_callback,alive_callback,visualization_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
# summary_callback() # print the timer summary
