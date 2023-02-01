using Plots
using Trixi
using Printf
using OrdinaryDiffEq



equations = Gaburro2D(1.0, 2.25*10^9, 1000, 9.81)

function initial_condition_test(x, t, equations::Gaburro2D)

    # # liquid domain
    # if((x[1]^2 + x[2]^2) <= 1)
    #     rho = 1000.0
    #     alpha = 1.0 - 10^(-3)
    #     v1 = -100.0 * x[1]
    #     v2 = 100.0 * x[2]
    # else
    #     rho = 1000.0
    #     v1 = 0.0
    #     v2 = 0.0
    #     alpha = 10^(-3)
    # end

  alpha_min = 10^(-1)

  rho = 1200.0
  v1 = 0.1
  v2 = 0.1
  alpha = alpha_min + (1.0 - alpha_min) * exp(-(x[1]^2 + x[2]^2)/1.0^2)
    
    return prim2cons(SVector(rho, v1, v2, alpha), equations)
end
  
initial_condition = initial_condition_test

# solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

surface_flux = (flux_lax_friedrichs, flux_nonconservative_gaburro)
volume_flux = (flux_central, flux_nonconservative_gaburro)

# solver = DGSEM(polydeg=1,
#   surface_flux=surface_flux, 
#   volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

@inline function gaburro_alpha(u, equations::Gaburro2D)
 return u[4]
end

basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                          alpha_max=1.0,
                                          alpha_min=0.001,
                                          alpha_smooth=true,
                                          variable=gaburro_alpha)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                  volume_flux_dg=volume_flux,
                                                  volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)



coordinates_min = (-3.0, -3.0)
coordinates_max = ( 3.0,  3.0)

cells_per_dimension = 64

mesh = StructuredMesh((cells_per_dimension,cells_per_dimension), coordinates_min, coordinates_max)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (0.0, 0.0076)
# tspan = (0.0, 0.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100

stepsize_callback = StepsizeCallback(cfl=0.8)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

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

  plotvar = plot_data["alpha"]
  # plotvar = plot_data["alpha_rho"]
  # rho, rhou, rhov, alpha = StructArrays.components(sol.u[end])

  # foo = rho # .* alpha

  # plotvar = ScalarPlotData2D(foo, semi)
  # plotvar = ScalarPlotData2D(x .* y, semi)

  # sol = ScalarPlotData2D(x .* y, semi))

  Plots.plot(plotvar,
    # clim=(0.0,1_500),
    colorbar_title="\ndensity",
    title=title,titlefontsize=10,
    dpi=300,
  )
  Plots.plot!(getmesh(plot_data),linewidth=0.5)

  mkpath("out")
  filename = joinpath("out", @sprintf("solution_%06d.png", timestep))
  Plots.savefig(filename)
end

function cons2myvar(u, equations)
  # assuming you are using compressible Euler
  rho = u[1]
  v1 = u[2]/u[1]
  v2 = u[3]/u[1]
  alpha = u[4]

  return SVector(rho * alpha)
end
Trixi.varnames(::typeof(cons2myvar), equations) = ("alpha_rho",) # , required to make it a tuple

visualization_callback = VisualizationCallback(
  plot_creator=my_save_plot,interval=100, clims=(0,1.1), show_mesh=true) #, solution_variables=cons2myvar)

callbacks = CallbackSet(stepsize_callback,alive_callback,visualization_callback)
# callbacks = CallbackSet(stepsize_callback,alive_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);

# summary_callback() # print the timer summary
