using Printf

"""
    T8codeMesh{NDIMS} <: AbstractMesh{NDIMS}

An unstructured curved mesh based on trees that uses the C library 't8code'
to manage trees and mesh refinement.
"""
mutable struct T8codeMesh{NDIMS, RealT<:Real, IsParallel, NDIMSP2, NNODES} <: AbstractMesh{NDIMS}
  cmesh                 :: Ptr{Cvoid} # cpointer to coarse mesh
  scheme                :: Ptr{Cvoid} # cpointer to element scheme
  forest                :: Ptr{Cvoid} # cpointer to forest
  is_parallel           :: IsParallel
  # ghost                 :: Ghost # Either Ptr{t8code_ghost_t} or Ptr{t8code_hex_ghost_t}
  # Coordinates at the nodes specified by the tensor product of 'nodes' (NDIMS times).
  # This specifies the geometry interpolation for each tree.
  tree_node_coordinates :: Array{RealT, NDIMSP2} # [dimension, i, j, k, tree]
  nodes                 :: SVector{NNODES, RealT}
  boundary_names        :: Array{Symbol, 2}      # [face direction, tree]
  current_filename      :: String
  unsaved_changes       :: Bool

  ncells                :: Int
  ninterfaces           :: Int
  nmortars              :: Int
  nboundaries           :: Int

  function T8codeMesh{NDIMS}(cmesh, scheme, forest, tree_node_coordinates, nodes, boundary_names,
                            current_filename, unsaved_changes) where NDIMS
    # if NDIMS == 2
    #   @assert t8code isa Ptr{t8code_t}
    # elseif NDIMS == 3
    #   @assert t8code isa Ptr{t8code_hex_t}
    # end

    # if mpi_isparallel()
    #   if !T8code.uses_mpi()
    #     error("t8code library does not support MPI")
    #   end
    #   is_parallel = Val(true)
    # else
    #   is_parallel = Val(false)
    # end
    is_parallel = Val(false)

    # ghost = ghost_new_t8code(t8code)

    mesh = new{NDIMS, eltype(tree_node_coordinates), typeof(is_parallel), NDIMS+2, length(nodes)}(
      cmesh, scheme, forest, is_parallel, tree_node_coordinates, nodes, boundary_names, current_filename, unsaved_changes)

    # Destroy 't8code' structs when the mesh is garbage collected.
    finalizer(function (mesh::T8codeMesh{NDIMS})
      trixi_t8_unref_forest(mesh.forest)
      finalize_t8code()
    end, mesh)

    return mesh
  end
end

SerialT8codeMesh{NDIMS} = T8codeMesh{NDIMS, <:Real, <:Val{false}}
# @inline mpi_parallel(mesh::SerialT8codeMesh) = Val(false)

@inline Base.ndims(::T8codeMesh{NDIMS}) where NDIMS = NDIMS
@inline Base.real(::T8codeMesh{NDIMS, RealT}) where {NDIMS, RealT} = RealT

# TODO: What should be returned in case of parallel processes? Local vs global.
@inline ntrees(mesh::T8codeMesh) = Int(t8_forest_get_num_global_trees(mesh.forest))
# @inline ncells(mesh::T8codeMesh) = Int(t8_forest_get_global_num_elements(mesh.forest))
@inline ncells(mesh::T8codeMesh) = Int(t8_forest_get_local_num_elements(mesh.forest))
@inline ninterfaces(mesh::T8codeMesh) = mesh.ninterfaces
@inline nmortars(mesh::T8codeMesh) = mesh.nmortars
@inline nboundaries(mesh::T8codeMesh) = mesh.nboundaries

function Base.show(io::IO, mesh::T8codeMesh)
  print(io, "T8codeMesh{", ndims(mesh), ", ", real(mesh), "}")
end

function Base.show(io::IO, :: MIME"text/plain", mesh::T8codeMesh)
  if get(io, :compact, false)
    show(io, mesh)
  else
    setup = [
             "#trees" => ntrees(mesh),
             "current #cells" => ncells(mesh),
             "polydeg" => length(mesh.nodes) - 1,
            ]
    summary_box(io, "T8codeMesh{" * string(ndims(mesh)) * ", " * string(real(mesh)) * "}", setup)
  end
end

"""
    T8codeMesh(trees_per_dimension; polydeg,
              mapping=nothing, faces=nothing, coordinates_min=nothing, coordinates_max=nothing,
              RealT=Float64, initial_refinement_level=0, periodicity=true, unsaved_changes=true)

Create a structured curved 'T8codeMesh' of the specified size.

There are three ways to map the mesh to the physical domain.
1. Define a 'mapping' that maps the hypercube '[-1, 1]^n'.
2. Specify a 'Tuple' 'faces' of functions that parametrize each face.
3. Create a rectangular mesh by specifying 'coordinates_min' and 'coordinates_max'.

Non-periodic boundaries will be called ':x_neg', ':x_pos', ':y_neg', ':y_pos', ':z_neg', ':z_pos'.

# Arguments
- 'trees_per_dimension::NTupleE{NDIMS, Int}': the number of trees in each dimension.
- 'polydeg::Integer': polynomial degree used to store the geometry of the mesh.
                      The mapping will be approximated by an interpolation polynomial
                      of the specified degree for each tree.
- 'mapping': a function of 'NDIMS' variables to describe the mapping that transforms
             the reference mesh ('[-1, 1]^n') to the physical domain.
- 'RealT::Type': the type that should be used for coordinates.
- 'initial_refinement_level::Integer': refine the mesh uniformly to this level before the simulation starts.
- 'periodicity': either a 'Bool' deciding if all of the boundaries are periodic or an 'NTuple{NDIMS, Bool}'
                 deciding for each dimension if the boundaries in this dimension are periodic.
- 'unsaved_changes::Bool': if set to 'true', the mesh will be saved to a mesh file.
"""
function T8codeMesh(trees_per_dimension; polydeg,
                   mapping, RealT=Float64, initial_refinement_level=0, periodicity=true, unsaved_changes=true)

  NDIMS = length(trees_per_dimension)

  # Convert periodicity to a Tuple of a Bool for every dimension
  if all(periodicity)
    # Also catches case where periodicity = true
    periodicity = ntuple(_->true, NDIMS)
  elseif !any(periodicity)
    # Also catches case where periodicity = false
    periodicity = ntuple(_->false, NDIMS)
  else
    # Default case if periodicity is an iterable
    periodicity = Tuple(periodicity)
  end

  conn = p4est_connectivity_new_brick(trees_per_dimension..., periodicity...)
  do_partition = 0
  cmesh = t8_cmesh_new_from_p4est(conn,t8_mpi_comm(),do_partition)
  p4est_connectivity_destroy(conn)

  scheme = t8_scheme_new_default_cxx()
  forest = t8_forest_new_uniform(cmesh,scheme,initial_refinement_level,0,mpi_comm().val)

  basis = LobattoLegendreBasis(RealT, polydeg)
  nodes = basis.nodes

  tree_node_coordinates = Array{RealT, NDIMS+2}(undef, NDIMS,
                                                ntuple(_ -> length(nodes), NDIMS)...,
                                                prod(trees_per_dimension))

  # Get cell length in reference mesh: Omega_ref = [-1,1]^2.
  dx = 2 / trees_per_dimension[1]
  dy = 2 / trees_per_dimension[2]

  num_local_trees = t8_cmesh_get_num_local_trees(cmesh)

  # Non-periodic boundaries
  boundary_names = fill(Symbol("---"), 2 * NDIMS, prod(trees_per_dimension))

  for itree = 1:num_local_trees
    veptr = t8_cmesh_get_tree_vertices(cmesh, itree-1)
    verts = unsafe_wrap(Array,veptr,(3,1 << NDIMS))

    # Calculate node coordinates of reference mesh.
    cell_x_offset = (verts[1,1] - 1/2*(trees_per_dimension[1]-1)) * dx
    cell_y_offset = (verts[2,1] - 1/2*(trees_per_dimension[2]-1)) * dy

    # u = [(mapping(cell_x_offset + dx * 1/2, cell_y_offset - dy * 1/2) - mapping(cell_x_offset - dx * 1/2, cell_y_offset - dy * 1/2))...,0.0]
    # v = [(mapping(cell_x_offset - dx * 1/2, cell_y_offset + dy * 1/2) - mapping(cell_x_offset - dx * 1/2, cell_y_offset - dy * 1/2))...,0.0]
    # w = [0.0,0.0,1.0]

    # vol = dot(cross(u,v),w)

    # if vol < 0
    #   sign = -1.0
    # else
    #   sign =  1.0
    # end

    for j in eachindex(nodes), i in eachindex(nodes)
      tree_node_coordinates[:, i, j, itree] .= mapping(cell_x_offset + dx * nodes[i]/2,
                                                       cell_y_offset + dy * nodes[j]/2)
    end

    if !periodicity[1]
      boundary_names[1, itree] = :x_neg
      boundary_names[2, itree] = :x_pos
    end

    if !periodicity[2]
      boundary_names[3, itree] = :y_neg
      boundary_names[4, itree] = :y_pos
    end
  end

  return T8codeMesh{NDIMS}(cmesh, scheme,forest, tree_node_coordinates, nodes, boundary_names, "", unsaved_changes)

end

function split_filename(filename)

  i = findlast(==('.'),filename)

  if isnothing(i)
    return filename
  end

  return filename[1:i-1],filename[i+1:end]

end

"""
    T8codeMesh{NDIMS}(meshfile::String;
                     mapping=nothing, polydeg=1, RealT=Float64,
                     initial_refinement_level=0, unsaved_changes=true)

Main mesh constructor for the `T8codeMesh` that imports an unstructured, conforming
mesh from an Abaqus mesh file (`.inp`). Each element of the conforming mesh parsed
from the `meshfile` is created as a [`p4est`](https://github.com/cburstedde/p4est)
tree datatype.

To create a curved unstructured mesh `T8codeMesh` two strategies are available:

- `p4est_mesh_from_hohqmesh_abaqus`: High-order, curved boundary information created by
                                     [`HOHQMesh.jl`](https://github.com/trixi-framework/HOHQMesh.jl) is
                                     available in the `meshfile`. The mesh polynomial degree `polydeg`
                                     of the boundaries is provided from the `meshfile`. The computation of
                                     the mapped tree coordinates is done with transfinite interpolation
                                     with linear blending similar to [`UnstructuredMesh2D`](@ref). Boundary name
                                     information is also parsed from the `meshfile` such that different boundary
                                     conditions can be set at each named boundary on a given tree.
- `p4est_mesh_from_standard_abaqus`: By default, with `mapping=nothing` and `polydeg=1`, this creates a
                                     straight-sided from the information parsed from the `meshfile`. If a mapping
                                     function is specified then it computes the mapped tree coordinates via polynomial
                                     interpolants with degree `polydeg`. The mesh created by this function will only
                                     have one boundary `:all`, as distinguishing different physical boundaries is
                                     non-trivial.

Note that the `mapping` and `polydeg` keyword arguments are only used by the `p4est_mesh_from_standard_abaqus`
function. The `p4est_mesh_from_hohqmesh_abaqus` function obtains the mesh `polydeg` directly from the `meshfile`
and constructs the transfinite mapping internally.

The particular strategy is selected according to the header present in the `meshfile` where
the constructor checks whether or not the `meshfile` was created with
[HOHQMesh.jl](https://github.com/trixi-framework/HOHQMesh.jl).
If the Abaqus file header is not present in the `meshfile` then the `P4estMesh` is created
with the function `p4est_mesh_from_standard_abaqus`.

The default keyword argument `initial_refinement_level=0` corresponds to a forest
where the number of trees is the same as the number of elements in the original `meshfile`.
Increasing the `initial_refinement_level` allows one to uniformly refine the base mesh given
in the `meshfile` to create a forest with more trees before the simulation begins.
For example, if a two-dimensional base mesh contains 25 elements then setting
`initial_refinement_level=1` creates an initial forest of `2^2 * 25 = 100` trees.

# Arguments
- `meshfile::String`: an uncurved Abaqus mesh file that can be imported by `p4est`.
- `mapping`: a function of `NDIMS` variables to describe the mapping that transforms
             the imported mesh to the physical domain. Use `nothing` for the identity map.
- `polydeg::Integer`: polynomial degree used to store the geometry of the mesh.
                      The mapping will be approximated by an interpolation polynomial
                      of the specified degree for each tree.
                      The default of `1` creates an uncurved geometry. Use a higher value if the mapping
                      will curve the imported uncurved mesh.
- `RealT::Type`: the type that should be used for coordinates.
- `initial_refinement_level::Integer`: refine the mesh uniformly to this level before the simulation starts.
- `unsaved_changes::Bool`: if set to `true`, the mesh will be saved to a mesh file.
"""
function T8codeMesh{NDIMS}(meshfile::String;
                          mapping=nothing, polydeg=1, RealT=Float64,
                          initial_refinement_level=0, unsaved_changes=true) where NDIMS
  # Prevent `t8code` from crashing Julia if the file doesn't exist.
  @assert isfile(meshfile)
  
  coordinates_min = (-1.0, -1.0)
  coordinates_max = ( 1.0,  1.0)

  # mapping = coordinates2mapping(coordinates_min, coordinates_max)

  # Create the mesh connectivity using `p4est`.
  # conn = p4est_connectivity_read_inp(meshfile)

  # do_partition = 0
  # cmesh = t8_cmesh_new_from_p4est(conn,t8_mpi_comm(),do_partition)
  # p4est_connectivity_destroy(conn)

  meshfile_prefix, meshfile_suffix = split_filename(meshfile)

  cmesh = t8_cmesh_from_msh_file(meshfile_prefix, 0, t8_mpi_comm(), 2, 0, 0)

  # error("debug")

  scheme = t8_scheme_new_default_cxx()
  forest = t8_forest_new_uniform(cmesh,scheme,initial_refinement_level,0,mpi_comm().val)

  basis = LobattoLegendreBasis(RealT, polydeg)
  nodes = basis.nodes

  # Get cell length in reference mesh: Omega_ref = [-1,1]^2.
  # dx = 2 / trees_per_dimension[1]
  # dy = 2 / trees_per_dimension[2]

  num_local_trees = t8_cmesh_get_num_local_trees(cmesh)

  tree_node_coordinates = Array{RealT, NDIMS+2}(undef, NDIMS,
                                                ntuple(_ -> length(nodes), NDIMS)...,
                                                num_local_trees)

  nodes_in = [-1.0, 1.0]
  matrix = polynomial_interpolation_matrix(nodes_in, nodes)
  data_in = Array{RealT, 3}(undef, 2, 2, 2)
  tmp1 = zeros(RealT, 2, length(nodes), length(nodes_in))

  for itree in 0:num_local_trees-1

    veptr = t8_cmesh_get_tree_vertices(cmesh, itree)
    verts = unsafe_wrap(Array,veptr,(3,1 << NDIMS))

    u = verts[:,2] - verts[:,1]
    v = verts[:,3] - verts[:,1]
    w = [0.0,0.0,1.0]

    vol = dot(cross(u,v),w)

    if vol < 0.0
      @warn "Discovered negative volumes in `cmesh`: vol = $vol"
    end

    # Tree vertices are stored what-not(?) order.
    @views data_in[:, 1, 1] .= verts[1:2,1]
    @views data_in[:, 2, 1] .= verts[1:2,2]
    @views data_in[:, 1, 2] .= verts[1:2,3]
    @views data_in[:, 2, 2] .= verts[1:2,4]

    # Interpolate corner coordinates to specified nodes.
    multiply_dimensionwise!(
      view(tree_node_coordinates, :, :, :, itree+1),
      matrix, matrix,
      data_in,
      tmp1
    )

  end

  map_node_coordinates!(tree_node_coordinates, mapping)

  # There's no simple and generic way to distinguish boundaries. Name all of them :all.
  boundary_names = fill(:all, 2 * NDIMS, num_local_trees)

  return T8codeMesh{NDIMS}(cmesh, scheme, forest, tree_node_coordinates, nodes, boundary_names, "", unsaved_changes)
end

function balance!(mesh::T8codeMesh{2}, init_fn=C_NULL)
  return nothing
end

function partition!(mesh::T8codeMesh{2}; allow_coarsening=true, weight_fn=C_NULL)
  return nothing
end

# Coarsen or refine marked cells and rebalance forest. Return a difference between
# old and new mesh.
function adapt!(mesh::T8codeMesh, indicators)

  old_levels = trixi_t8_get_local_element_levels(mesh.forest)

  forest_cached = trixi_t8_adapt_new(mesh.forest, indicators)

  new_levels = trixi_t8_get_local_element_levels(forest_cached)

  differences = trixi_t8_get_difference(old_levels, new_levels)

  mesh.forest = forest_cached

  return differences
end
