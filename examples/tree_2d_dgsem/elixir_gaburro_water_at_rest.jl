using OrdinaryDiffEq
using Trixi
using Printf


equations = Gaburro2D(1.0, 2.78*10^5, 1000.0, 9.81)

function initial_condition_rest(x, t, equations::Gaburro2D)
    # liquid domain
    rho = equations.rho_0 * exp(-(equations.gravity * equations.rho_0/equations.k0) *(x[2] - 1.0))
    v1 = 0.0
    v2 = 0.0
    alpha = 1.0
    
    return prim2cons(SVector(rho, v1, v2, alpha), equations)
end
  
initial_condition = initial_condition_rest

function boundary_condition_wallG(u_inner, orientation, direction, x, t, surface_flux_function,
    equations::Gaburro2D)
  
    if direction in (1, 2) # x direction
      u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])
    else # y direction
      u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4])
    end
  
    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
      flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
      flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end
  
    return flux
  end
  

boundary_conditions = (x_neg=boundary_condition_wallG,
                       x_pos=boundary_condition_wallG,
                       y_neg=boundary_condition_wallG,
                       y_pos=boundary_condition_wallG,)
  
volume_flux = (flux_central, flux_nonconservative_gaburro)
surface_flux=(flux_lax_friedrichs, flux_nonconservative_gaburro)
#solver = DGSEM(polydeg=3, surface_flux=(flux_lax_friedrichs, flux_nonconservative_gaburro),
#                 volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                          alpha_max=0.5,
                                          alpha_min=0.001,
                                          alpha_smooth=true,
                                          variable=density)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                  volume_flux_dg=volume_flux,
                                                  volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-0.5, 0.0) # minimum coordinates (min(x), min(y))
coordinates_max = ( 0.5, 1.0) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh with periodic boundaries
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=30_000, periodicity=(false,false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, 
                source_terms=source_terms_gravity, boundary_conditions=boundary_conditions)

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100

alive_callback = AliveCallback(analysis_interval=analysis_interval)

stepsize_callback = StepsizeCallback(cfl=0.8)

function save_my_plot(plot_data, variable_names;
  show_mesh=false, plot_arguments=Dict{Symbol,Any}(),
  time=nothing, timestep=nothing)

  # Gather subplots
  plots = []
  for v in variable_names
    push!(plots, Plots.plot(plot_data[v]; plot_arguments...))
  end
  if show_mesh
    push!(plots, Plots.plot(getmesh(plot_data); plot_arguments...))
  end

  pressure_matrix = equations.k0 .* plot_data.data[1]
  pressure_matrix = pressure_matrix .- equations.k0
  push!(plots, Plots.plot(heatmap(plot_data.x, plot_data.y, pressure_matrix), title = "pressure", width=10, height=10))

  # Create plot
  Plots.plot(plots...,)

  # Determine filename and save plot
  filename = joinpath("out", @sprintf("solution_%06d.png", timestep))
  Plots.savefig(filename)
end

visualization_callback = VisualizationCallback(; interval=1000,
                          solution_variables=cons2prim,
                          #variable_names=["rho"],
                          show_mesh=false,
                          plot_data_creator=PlotData2D,
                          plot_creator=save_my_plot,
                          )

callbacks = CallbackSet(stepsize_callback, visualization_callback, alive_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary