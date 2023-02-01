using OrdinaryDiffEq
using Revise
using Trixi
using Printf
using Plots

equations = Gaburro2D(1.0, 2.25*10^9, 1000.0, 0.0)

function initial_condition_test(x, t, equations::Gaburro2D)
  
    # liquid domain
    if((x[1]^2 + x[2]^2) <= 1)
        rho = 1000.0
        alpha = 1.0 - 10^(-3)
        v1 = -100.0 * x[1]
        v2 = 100.0 * x[2]
    else
        rho = 1000.0
        v1 = 0.0
        v2 = 0.0
        alpha = 10^(-3)
    end
    
    return prim2cons(SVector(rho, v1, v2, alpha), equations)
end

initial_condition = initial_condition_test

@inline function gaburro_alpha_rho(u, equations::Gaburro2D)
 return u[1]
end

volume_flux = (flux_central, flux_nonconservative_gaburro)
surface_flux=(flux_lax_friedrichs, flux_nonconservative_gaburro)


# basis = LobattoLegendreBasis(3)
# indicator_sc = IndicatorHennemannGassner(equations, basis,
#                                           alpha_max=1.0,
#                                           alpha_min=0.001,
#                                           alpha_smooth=true,
#                                           variable=gaburro_alpha_rho)
# volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
#                                                   volume_flux_dg=volume_flux,
#                                                   volume_flux_fv=surface_flux)

volume_integral=VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)


coordinates_min = (-3.0, -3.0) # minimum coordinates (min(x), min(y))
coordinates_max = ( 3.0,  3.0) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh with periodic boundaries
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=8,
                n_cells_max=30_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (0.0, 0.0076)
# tspan = (0.0, 0.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100

alive_callback = AliveCallback(analysis_interval=analysis_interval)

stepsize_callback = StepsizeCallback(cfl=0.5)

# function save_my_plot(plot_data, variable_names;
#     show_mesh=false, plot_arguments=Dict{Symbol,Any}(),
#     time=nothing, timestep=nothing)
# 
#   # Gather subplots
#   plots = []
#   for v in variable_names
#     push!(plots, Plots.plot(plot_data[v]; plot_arguments...))
#   end
#   if show_mesh
#     push!(plots, Plots.plot(getmesh(plot_data); plot_arguments...))
#   end
# 
#   pressure_matrix = equations.k0 .* plot_data.data[1]
#   pressure_matrix = pressure_matrix .- equations.k0
#   push!(plots, Plots.plot(heatmap(plot_data.x, plot_data.y, pressure_matrix), title = "pressure", width=10, height=10))
# 
#   # Create plot
#   Plots.plot(plots...,)
# 
#   # Determine filename and save plot
#   filename = joinpath("out", @sprintf("solution_%06d.png", timestep))
#   Plots.savefig(filename)
# end

function my_save_plot(plot_data, variable_names;
                   show_mesh=true, plot_arguments=Dict{Symbol,Any}(),
                   time=nothing, timestep=nothing)

  title = @sprintf("2D KHI | Trixi.jl | 2nd-order DG | Gaburro: t = %3.2f", time)

  # println("plot_data")
  # println(plot_data["rho"].plot_data.x)

  # x = plot_data["rho"].plot_data.x
  # y = plot_data["rho"].plot_data.y

  # rho = plot_data["rho"].plot_data.data
  # alpha = plot_data["alpha"].plot_data.data

  # sol = plot_data["rho"] # .* plot_data["alpha"]

  # sol = plot_data["rho"].plot_data

  # plotvar = plot_data["alpha"]
  # plotvar = plot_data["rho"]
  plotvar = plot_data["alpha_rho"]
  # plotvar = plot_data["alpha_rho"]
  # rho, rhou, rhov, alpha = StructArrays.components(sol.u[end])

  # foo = rho # .* alpha

  # plotvar = ScalarPlotData2D(foo, semi)
  # plotvar = ScalarPlotData2D(x .* y, semi)

  # sol = ScalarPlotData2D(x .* y, semi))

  Plots.plot(plotvar,
    clim=(0.0,1200.0),
    colorbar_title="\ndensity",
    title=title,titlefontsize=10,
    dpi=300,
  )

  # Plots.plot!(getmesh(plot_data),linewidth=0.5)

  mkpath("out")
  filename = joinpath("out", @sprintf("solution_%06d.png", timestep))
  Plots.savefig(filename)
end

# Convert conservative variables to primitive
function cons2cons(u, equations::Gaburro2D)
  alpha_rho, alpha_rho_v1, alpha_rho_v2, alpha = u

  # rho = alpha_rho/alpha
  # v1 = alpha_rho_v1 / alpha_rho
  # v2 = alpha_rho_v2 / alpha_rho
    
  return SVector(alpha_rho, alpha_rho_v1, alpha_rho_v2, alpha)
end

Trixi.varnames(::typeof(cons2cons), equations) = ("alpha_rho", "alpha_rho_v1", "alpha_rho_v2", "alpha")
 

#visualization_callback = VisualizationCallback(plot_creator=my_save_plot,interval=10, clims=(0,1.1), show_mesh=true)
visualization_callback = VisualizationCallback(; interval=500,
                            # solution_variables=cons2prim,
                            solution_variables=cons2cons,
                            #variable_names=["rho"],
                            show_mesh=false,
                            # plot_data_creator=PlotData2D,
                            plot_creator=my_save_plot,
                            #plot_creator=save_my_plot,
                            )

callbacks = CallbackSet(stepsize_callback, visualization_callback, alive_callback)

###############################################################################
# run the simulation

#stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-6, 5.0e-6),
#                                                     variables=(Trixi.density,density))

sol = solve(ode, CarpenterKennedy2N54( williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

# deg = 3
# level = 5
# sol_all_var = sol[2]
# sol_all_var_0 = sol[1]
# 
# rho_vec = zeros((deg + 1)^2 * 2^level * 2^level)
# rho0_vec = zeros((deg + 1)^2 * 2^level * 2^level)
# 
# for i = 1:((deg + 1)^2 * 2^level * 2^level)
#     rho_vec[i] = sol_all_var[(i-1)*4 + 1]/sol_all_var[(i-1)*4 + 4]
#     rho0_vec[i] = sol_all_var_0[(i-1)*4 + 1]/sol_all_var_0[(i-1)*4 + 4]
# end
# 
# x_LGL, w_LGL = Trixi.gauss_lobatto_nodes_weights(deg+1)
# 
# integral_rho = 0
# integral_rho0 = 0
# for i = 1:(deg + 1)^2:length(rho0_vec)
#     for k in 1:(deg + 1)
#         for l in 1:(deg + 1)
#             global integral_rho += rho_vec[i + (k-1) + (l-1)] * w_LGL[l] * w_LGL[k] #* (coordinates_max-coordinates_min)/(2^level)*0.5
#             global integral_rho0 += rho0_vec[i + (k-1) + (l-1)] * w_LGL[l] * w_LGL[k] #* (coordinates_max-coordinates_min)/(2^level)*0.5
#         end
#     end
# end
# 
# println("Fehler Integral über erste Variable:  ", abs(integral_rho - integral_rho0))
