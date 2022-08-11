"""
    T8codeMesh{NDIMS} <: AbstractMesh{NDIMS}

An unstructured curved mesh based on trees that uses the C library 't8code'
to manage trees and mesh refinement.
"""
mutable struct T8codeMesh{NDIMS, RealT<:Real, IsParallel, NDIMSP2, NNODES} <: AbstractMesh{NDIMS}
  cmesh                 :: Ptr{Cvoid} # cpointer to coarse mesh
  scheme                :: Ptr{Cvoid} # cpointer to element scheme
  forest                :: Ptr{Cvoid} # cpointer to forest
  forest_cached         :: Ptr{Cvoid} # cpointer to a cached forest
  is_parallel           :: IsParallel
  # ghost                 :: Ghost # Either Ptr{t8code_ghost_t} or Ptr{t8code_hex_ghost_t}
  # Coordinates at the nodes specified by the tensor product of 'nodes' (NDIMS times).
  # This specifies the geometry interpolation for each tree.
  tree_node_coordinates :: Array{RealT, NDIMSP2} # [dimension, i, j, k, tree]
  nodes                 :: SVector{NNODES, RealT}
  boundary_names        :: Array{Symbol, 2}      # [face direction, tree]
  current_filename      :: String
  unsaved_changes       :: Bool

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
      cmesh, scheme, forest, C_NULL, is_parallel, tree_node_coordinates, nodes, boundary_names, current_filename, unsaved_changes)

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
@inline ncells(mesh::T8codeMesh) = Int(t8_forest_get_global_num_elements(mesh.forest))
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
             Use only one of 'mapping', 'faces' and 'coordinates_min'/'coordinates_max'.
- 'faces::NTuple{2*NDIMS}': a tuple of '2 * NDIMS' functions that describe the faces of the domain.
                            Each function must take 'NDIMS-1' arguments.
                            'faces[1]' describes the face onto which the face in negative x-direction
                            of the unit hypercube is mapped. The face in positive x-direction of
                            the unit hypercube will be mapped onto the face described by 'faces[2]'.
                            'faces[3:4]' describe the faces in positive and negative y-direction respectively
                            (in 2D and 3D).
                            'faces[5:6]' describe the faces in positive and negative z-direction respectively (in 3D).
                            Use only one of 'mapping', 'faces' and 'coordinates_min'/'coordinates_max'.
- 'coordinates_min': vector or tuple of the coordinates of the corner in the negative direction of each dimension
                     to create a rectangular mesh.
                     Use only one of 'mapping', 'faces' and 'coordinates_min'/'coordinates_max'.
- 'coordinates_max': vector or tuple of the coordinates of the corner in the positive direction of each dimension
                     to create a rectangular mesh.
                     Use only one of 'mapping', 'faces' and 'coordinates_min'/'coordinates_max'.
- 'RealT::Type': the type that should be used for coordinates.
- 'initial_refinement_level::Integer': refine the mesh uniformly to this level before the simulation starts.
- 'periodicity': either a 'Bool' deciding if all of the boundaries are periodic or an 'NTuple{NDIMS, Bool}'
                 deciding for each dimension if the boundaries in this dimension are periodic.
- 'unsaved_changes::Bool': if set to 'true', the mesh will be saved to a mesh file.
"""
function T8codeMesh(trees_per_dimension; polydeg,
                   mapping=nothing, faces=nothing, coordinates_min=nothing, coordinates_max=nothing,
                   RealT=Float64, initial_refinement_level=0, periodicity=true, unsaved_changes=true)

  @assert (
    (coordinates_min === nothing) === (coordinates_max === nothing)
  ) "Either both or none of coordinates_min and coordinates_max must be specified"

  @assert count(i -> i !== nothing,
    (mapping, faces, coordinates_min)
  ) == 1 "Exactly one of mapping, faces and coordinates_min/max must be specified"

  # Extract mapping
  if faces !== nothing
    validate_faces(faces)
    mapping = transfinite_mapping(faces)
  elseif coordinates_min !== nothing
    mapping = coordinates2mapping(coordinates_min, coordinates_max)
  end

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

  basis = LobattoLegendreBasis(RealT, polydeg)
  nodes = basis.nodes
  tree_node_coordinates = Array{RealT, NDIMS+2}(undef, NDIMS,
                                                ntuple(_ -> length(nodes), NDIMS)...,
                                                prod(trees_per_dimension))
  calc_tree_node_coordinates!(tree_node_coordinates, nodes, mapping, trees_per_dimension)

  # p4est_connectivity_new_brick has trees in Z-order, so use our own function for this
  # conn = connectivity_structured(trees_per_dimension..., periodicity)

  # println("")
  # println("conn   = ", conn)

  # cmesh = t8_cmesh_new_from_p4est(conn,t8_mpi_comm(),0)

  scheme = t8_scheme_new_default_cxx()
  cmesh = t8_cmesh_new_periodic(t8_mpi_comm(),NDIMS)
  forest = t8_forest_new_uniform(cmesh,scheme,initial_refinement_level,0,mpi_comm().val)

  # println("scheme = ", scheme)
  # println("cmesh  = ", cmesh)
  # println("forest = ", forest)
  # println("")

  # t8code = new_t8code(conn, initial_refinement_level)
  # cmesh,forest,scheme = trixi_t8code_new_uniform_mesh(conn, initial_refinement_level)

  # Non-periodic boundaries
  boundary_names = fill(Symbol("---"), 2 * NDIMS, prod(trees_per_dimension))

  structured_boundary_names!(boundary_names, trees_per_dimension, periodicity)

  return T8codeMesh{NDIMS}(cmesh, scheme,forest, tree_node_coordinates, nodes, boundary_names, "", unsaved_changes)
end

function balance!(mesh::T8codeMesh{2}, init_fn=C_NULL)
  # t8code_balance(mesh.t8code, T8CODE_DQUAD_CONNECT_FACE, init_fn)
  return nothing
end

function partition!(mesh::T8codeMesh{2}; allow_coarsening=true, weight_fn=C_NULL)
  # t8code_partition(mesh.t8code, Int(allow_coarsening), weight_fn)
  return nothing
end

# Coarsen or refine marked cells and rebalance forest. Return a list of all
# cells that have been adapted during coarsening/refinement or rebalancing.
function adapt!(mesh::T8codeMesh, indicators)

  println("## begin t8code_adapt!")
  old_levels = trixi_t8_get_local_element_levels(mesh.forest)

  mesh.forest_cached = trixi_t8_adapt_new(mesh.forest, indicators)

  new_levels = trixi_t8_get_local_element_levels(mesh.forest_cached)

  # differences = trixi_t8_get_difference(mesh.forest, mesh.forest_cached)
  differences = trixi_t8_get_difference(old_levels, new_levels)

  trixi_t8_unref_forest(mesh.forest)
  mesh.forest = mesh.forest_cached
  println("## end t8code_adapt!")

  return differences
end

# Coarsen marked cells if the forest will stay balanced.
# Return a list of all cells that have been coarsened.
# function coarsen!(mesh::T8codeMesh)
#   return nothing
# end
