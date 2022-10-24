"""
    T8codeHybridMesh{NDIMS} <: AbstractMesh{NDIMS}

An unstructured curved hybrid mesh based on trees that uses the C library 't8code'
to manage trees and mesh refinement.
"""
mutable struct T8codeHybridMesh{NDIMS, RealT<:Real, IsParallel} <: AbstractMesh{NDIMS}
  cmesh                 :: Ptr{Cvoid} # cpointer to coarse mesh
  scheme                :: Ptr{Cvoid} # cpointer to element scheme
  forest                :: Ptr{Cvoid} # cpointer to forest
  is_parallel           :: IsParallel

  # ghost                 :: Ghost # Either Ptr{t8code_ghost_t} or Ptr{t8code_hex_ghost_t}
  # Coordinates at the nodes specified by the tensor product of 'nodes' (NDIMS times).
  # This specifies the geometry interpolation for each tree.
  
  # tree_node_coordinates # :: Array{RealT, NDIMS+2} # [dimension, i, j, k, tree]

  # nodes                 # :: SVector{NNODES, RealT}

  boundary_names        :: Array{Symbol, 2}      # [face direction, tree]

  # current_filename      :: String
  # unsaved_changes       :: Bool

  element_types         :: Set{Symbol}
  mesh_datas            :: Array{MeshData, 1}

  ncells                :: Int
  ninterfaces           :: Int
  nmortars              :: Int
  nboundaries           :: Int

  function T8codeHybridMesh{NDIMS}(cmesh, scheme, forest) where NDIMS
    is_parallel = Val(false)

    mesh = new{NDIMS, Float64, typeof(is_parallel)}(cmesh, scheme, forest, is_parallel)

    # Destroy 't8code' structs when the mesh is garbage collected.
    finalizer(function (mesh::T8codeHybridMesh{NDIMS})
      trixi_t8_unref_forest(mesh.forest)
      finalize_t8code()
    end, mesh)

    return mesh
  end
end

function T8codeHybridMesh{NDIMS}(cmesh, scheme, forest, boundary_names) where NDIMS

  mesh = T8codeHybridMesh{NDIMS}(cmesh, scheme, forest)

  mesh.boundary_names = boundary_names

  # mesh.nodes = nodes
  # mesh.unsaved_changes = unsaved_changes
  # mesh.current_filename = current_filename
  # mesh.tree_node_coordinates = tree_node_coordinates

  return mesh
end

SerialT8codeHybridMesh{NDIMS} = T8codeHybridMesh{NDIMS, <:Val{false}}
# @inline mpi_parallel(mesh::SerialT8codeHybridMesh) = Val(false)

@inline Base.ndims(::T8codeHybridMesh{NDIMS}) where NDIMS = NDIMS
@inline Base.real(::T8codeHybridMesh{NDIMS, RealT}) where {NDIMS, RealT} = RealT

# TODO: What should be returned in case of parallel processes? Local vs global.
@inline ntrees(mesh::T8codeHybridMesh) = Int(t8_forest_get_num_global_trees(mesh.forest))
# @inline ncells(mesh::T8codeHybridMesh) = Int(t8_forest_get_global_num_elements(mesh.forest))
@inline ncells(mesh::T8codeHybridMesh) = Int(t8_forest_get_local_num_elements(mesh.forest))
@inline ninterfaces(mesh::T8codeHybridMesh) = mesh.ninterfaces
@inline nmortars(mesh::T8codeHybridMesh) = mesh.nmortars
@inline nboundaries(mesh::T8codeHybridMesh) = mesh.nboundaries

function Base.show(io::IO, mesh::T8codeHybridMesh)
  print(io, "T8codeHybridMesh{", ndims(mesh), ", ", real(mesh), "}")
end

function Base.show(io::IO, :: MIME"text/plain", mesh::T8codeHybridMesh)
  if get(io, :compact, false)
    show(io, mesh)
  else
    setup = [
             "#trees" => ntrees(mesh),
             "current #cells" => ncells(mesh),
             # "polydeg" => length(mesh.nodes) - 1,
            ]
    summary_box(io, "T8codeHybridMesh{" * string(ndims(mesh)) * ", " * string(real(mesh)) * "}", setup)
  end
end

"""
    T8codeHybridMesh(cmesh; polydeg,
              mapping=nothing,
              RealT=Float64, initial_refinement_level=0)

Create an un-structured curved 'T8codeHybridMesh' according to a given `cmesh`. 

# Arguments
- 'cmesh::Ptr{Cvoid}': C-pointer to a `cmesh` datastructure.
- 'polydeg::Integer': polynomial degree used to store the geometry of the mesh.
                      The mapping will be approximated by an interpolation polynomial
                      of the specified degree for each tree.
- 'mapping': a function of 'NDIMS' variables to describe the mapping that transforms
             the reference mesh ('[-1, 1]^n') to the physical domain.
- 'RealT::Type': the type that should be used for coordinates.
- 'initial_refinement_level::Integer': refine the mesh uniformly to this level before the simulation starts.
"""

function T8codeHybridMesh(cmesh :: Ptr{Cvoid}; NDIMS, polydeg, mapping, RealT=Float64, initial_refinement_level=0)

  scheme = t8_scheme_new_default_cxx()
  forest = t8_forest_new_uniform(cmesh,scheme,initial_refinement_level,0,mpi_comm().val)

  num_local_trees = t8_cmesh_get_num_local_trees(cmesh)

  # Building Gauss-Lobatto basis for a quadrangle.
  # r1D, w1D = StartUpDG.gauss_lobatto_quad(0, 0, polydeg)
  # rq, sq = vec.(StartUpDG.meshgrid(r1D))
  # wr, ws = vec.(StartUpDG.meshgrid(w1D))
  # wq = @. wr * ws
  # rd = RefElemData(Quad(), polydeg; quad_rule_vol = (rq, sq, wq), quad_rule_face = (r1D, w1D))
    
  # tree_node_coordinates = Array{RealT, NDIMS+2}(undef, NDIMS,
  #                                               ntuple(_ -> length(nodes), NDIMS)...,
  #                                               num_local_trees)

  # Non-periodic boundaries
  # boundary_names = fill(Symbol("---"), 2 * NDIMS, prod(trees_per_dimension))

  for itree = 1:num_local_trees

    eclass = t8_eclass_to_element_type[t8_cmesh_get_tree_class(cmesh, itree-1)]
    nverts = t8_eclass_num_vertices[eclass]
    tverts = unsafe_wrap(Array,t8_cmesh_get_tree_vertices(cmesh, itree-1),(3,nverts))

    println(eclass)
    display(tverts)
    println("")
    
    # mapping = coordinates2mapping((tverts[1,1],tverts[2,1]),(tverts[1,4],tverts[2,4]))

    # (VX,VY),EToV = StartUpDG.uniform_mesh(Quad(),1,1)
    # md = MeshData((VX, VY), EToV, rd)

    # @unpack x,y = md
    # for i in eachindex(x)
    #   x[i],y[i] = mapping(x[i],y[i])
    # end

    # md = MeshData(rd, md, x, y)

    # println(x)
    # println(y)
    # println("")

    # for j in eachindex(nodes), i in eachindex(nodes)
    #   tree_node_coordinates[:, i, j, itree] .= mapping(cell_x_offset + dx * nodes[i]/2,
    #                                                    cell_y_offset + dy * nodes[j]/2)
    # end

    # if !periodicity[1]
    #   boundary_names[1, itree] = :x_neg
    #   boundary_names[2, itree] = :x_pos
    # end

    # if !periodicity[2]
    #   boundary_names[3, itree] = :y_neg
    #   boundary_names[4, itree] = :y_pos
    # end
  end

  # error("debug")

  # There's no simple and generic way to distinguish boundaries. Name all of them :all.
  boundary_names = fill(:all, 2 * NDIMS, num_local_trees)

  nfaces = trixi_t8_count_faces(forest)

  println(nfaces)

  levels = Vector{Int}(undef, t8_forest_get_local_num_elements(forest))
  shapes = Vector{Symbol}(undef, t8_forest_get_local_num_elements(forest))

  mapP = Vector{Int}(undef, nfaces)
  mapM = Vector{Int}(undef, nfaces)
  orientation = Vector{Int}(undef, nfaces)

  trixi_t8_fill_mesh_info!(forest, shapes, levels, mapP, mapM, orientation)

  println(levels)
  println(mapP)
  println(mapM)
  println(orientation)

  # trixi_t8_fill_mesh_info(mesh.forest, elements, interfaces, mortars, boundaries, mesh.boundary_names)

  return T8codeHybridMesh{NDIMS}(cmesh, scheme, forest, boundary_names)

end

function balance!(mesh::T8codeHybridMesh{2}, init_fn=C_NULL)
  return nothing
end

function partition!(mesh::T8codeHybridMesh{2}; allow_coarsening=true, weight_fn=C_NULL)
  return nothing
end

# Coarsen or refine marked cells and rebalance forest. Return a difference between
# old and new mesh.
function adapt!(mesh::T8codeHybridMesh, indicators)

  # old_levels = trixi_t8_get_local_element_levels(mesh.forest)

  # forest_cached = trixi_t8_adapt_new(mesh.forest, indicators)

  # new_levels = trixi_t8_get_local_element_levels(forest_cached)

  # differences = trixi_t8_get_difference(old_levels, new_levels)

  # mesh.forest = forest_cached

  # return differences

end
