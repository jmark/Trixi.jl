@muladd begin

# This method is called when a SemidiscretizationHyperbolic is constructed.
# It constructs the basic `cache` used throughout the simulation to compute
# the RHS etc.
function create_cache(mesh::T8codeMesh, equations::AbstractEquations, dg::DG, ::Any, ::Type{uEltype}) where {uEltype<:Real}
  # Make sure to balance the 't8code' before creating any containers
  # in case someone has tampered with the 't8code' after creating the mesh.
  # balance!(mesh)

  _ = count_required_surfaces(mesh)

  elements   = init_elements(mesh, equations, dg.basis, uEltype)
  interfaces = init_interfaces(mesh, equations, dg.basis, elements)
  boundaries = init_boundaries(mesh, equations, dg.basis, elements)
  mortars    = init_mortars(mesh, equations, dg.basis, elements)

  _ = init_surfaces!(interfaces, mortars, boundaries, mesh)

  cache = (; elements, interfaces, boundaries, mortars)

  # Add specialized parts of the cache required to compute the volume integral etc.
  cache = (; cache..., create_cache(mesh, equations, dg.volume_integral, dg, uEltype)...)
  cache = (; cache..., create_cache(mesh, equations, dg.mortar, uEltype)...)

  return cache
end

include("containers.jl")
include("containers_2d.jl")
include("dg_2d.jl")

# Not implemented yet.
# include("dg_3d.jl")
# include("dg_parallel.jl")

end # @muladd
