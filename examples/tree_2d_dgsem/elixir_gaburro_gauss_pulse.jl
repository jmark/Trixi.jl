using OrdinaryDiffEq
using Revise
using Trixi
using Printf

equations = Gaburro2D(1.0, 2.25*10^2, 1.0, 9.81)


function initial_condition_gauss_pulse(x, t, equations::Gaburro2D)
    
    if(((x[1]^2 + x[2]^2) <= 2))
        # liquid domain   
        rho = 1.0 * exp(-1/(1.0)^2*(x[1]^2 + x[2]^2)) + 1.0
        v1 = 0.0
        v2 = 0.0
        alpha = 1.0
    else
        rho = 1.0 * exp(-1/(1.0)^2*(x[1]^2 + x[2]^2)) + 1.0
        v1 = 0.0
        v2 = 0.0
        alpha = 0.8
    end

    return prim2cons(SVector(rho, v1, v2, alpha), equations)
end

function initial_condition_blast_wave(x, t, equations::Gaburro2D)
    # Modified From Hennemann & Gassner JCP paper 2020 (Sec. 6.3) -> "medium blast wave"
    # Set up polar coordinates
    inicenter = SVector(0.0, 0.0)
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    phi = atan(y_norm, x_norm)
    sin_phi, cos_phi = sincos(phi)
  
    # Calculate primitive variables
    rho = r > 0.5 ? 1.0 : 1.1691
    v1  = r > 0.5 ? 0.0 : 0.1882 * cos_phi
    v2  = r > 0.5 ? 0.0 : 0.1882 * sin_phi
    alpha = r > 0.5 ? 0.8 : 1.0
  
    return prim2cons(SVector(rho, v1, v2, alpha), equations)
  end
  
initial_condition = initial_condition_gauss_pulse


volume_flux = (flux_central, flux_nonconservative_gaburro)
surface_flux=(flux_lax_friedrichs, flux_nonconservative_gaburro)
#solver = DGSEM(polydeg=1, surface_flux=surface_flux,
 #                volume_integral=VolumeIntegralFluxDifferencing(volume_flux))


@inline function gaburro_alpha_rho(u, equations::Gaburro2D)
 return u[1]
end

basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                          alpha_max=0.5,
                                          alpha_min=0.001,
                                          alpha_smooth=true,
                                          variable=gaburro_alpha_rho)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                  volume_flux_dg=volume_flux,
                                                  volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-2.0, -2.0) # minimum coordinates (min(x), min(y))
coordinates_max = ( 2.0, 2.0) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh with periodic boundaries
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=5,
                n_cells_max=10_000,
                periodicity=false)

#boundary_conditions = BoundaryConditionDirichlet(initial_condition_gauss_pulse)                
#boundary_conditions = boundary_condition_do_nothing
#boundary_conditions = boundary_condition_wall

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms = source_terms_gravity,
                                    boundary_conditions=boundary_condition_wall)

tspan = (0.0, 1.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100

alive_callback = AliveCallback(analysis_interval=analysis_interval)

stepsize_callback = StepsizeCallback(cfl=0.4)

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

visualization_callback = VisualizationCallback(; interval=100,
                            solution_variables=cons2prim,
                            #variable_names=["rho"],
                            show_mesh=false,
                            plot_data_creator=PlotData2D,
                            #plot_creator=save_my_plot,
                            )

callbacks = CallbackSet(alive_callback, visualization_callback, stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary
