@muladd begin

# Initialize data structures in element container.
function init_elements!(elements, mesh::T8codeMesh{2}, basis::LobattoLegendreBasis)
  @unpack node_coordinates, jacobian_matrix,
          contravariant_vectors, inverse_jacobian = elements

  calc_node_coordinates!(node_coordinates, mesh, basis)

  for element in 1:ncells(mesh)
    calc_jacobian_matrix!(jacobian_matrix, element, node_coordinates, basis)

    calc_contravariant_vectors!(contravariant_vectors, element, jacobian_matrix)

    calc_inverse_jacobian!(inverse_jacobian, element, jacobian_matrix)
  end

  return nothing
end

# Interpolate tree_node_coordinates to each quadrant at the nodes of the specified basis.
function calc_node_coordinates!(node_coordinates,
                                mesh::T8codeMesh{2},
                                basis::LobattoLegendreBasis)
  # Hanging nodes will cause holes in the mesh if its polydeg is higher
  # than the polydeg of the solver.
  @assert length(basis.nodes) >= length(mesh.nodes) "The solver can't have a lower polydeg than the mesh."

  calc_node_coordinates!(node_coordinates, mesh, basis.nodes)
end

# Interpolate tree_node_coordinates to each quadrant at the specified nodes.
function calc_node_coordinates!(node_coordinates,
                                mesh::T8codeMesh{2},
                                nodes::AbstractVector)
  # We use `StrideArray`s here since these buffers are used in performance-critical
  # places and the additional information passed to the compiler makes them faster
  # than native `Array`s.
  tmp1    = StrideArray(undef, real(mesh),
                        StaticInt(2), static_length(nodes), static_length(mesh.nodes))
  matrix1 = StrideArray(undef, real(mesh),
                        static_length(nodes), static_length(mesh.nodes))
  matrix2 = similar(matrix1)
  baryweights_in = barycentric_weights(mesh.nodes)

  num_local_trees = t8_forest_get_num_local_trees(mesh.forest)

  current_index = 0
  for itree = 0:num_local_trees-1

    tree_class = t8_forest_get_tree_class(mesh.forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(mesh.forest, tree_class);

    num_elements_in_tree = t8_forest_get_tree_num_elements(mesh.forest, itree);

    for ielement = 0:num_elements_in_tree-1
      current_index += 1

      # TODO: Make this more general.
      element = t8_forest_get_element_in_tree(mesh.forest, itree, ielement)
      element_level  = t8_element_level(eclass_scheme, element)
      element_length = t8code_quadrant_len(element_level) / t8code_root_len
      
      centroid = Array{Float64}(undef,3);      # /* Will hold the element midpoint. */
      t8_forest_element_centroid(mesh.forest, itree, element, pointer(centroid));

      nodes_out_x = 2*(element_length * 1/2 * nodes .+ centroid[1]) .- 1.0
      nodes_out_y = 2*(element_length * 1/2 * nodes .+ centroid[2]) .- 1.0
      # end TODO

      polynomial_interpolation_matrix!(matrix1, mesh.nodes, nodes_out_x, baryweights_in)
      polynomial_interpolation_matrix!(matrix2, mesh.nodes, nodes_out_y, baryweights_in)

      multiply_dimensionwise!(
        view(node_coordinates, :, :, :, ielement+1),
        matrix1, matrix2,
        view(mesh.tree_node_coordinates, :, :, :, itree+1),
        tmp1
      )
    end
  end

  return node_coordinates
end

# Initialize node_indices of interface container.
@inline function init_interface_node_indices!(interfaces::T8codeInterfaceContainer{2},
                                              faces, orientation, interface_id)
  # Iterate over primary and secondary element.
  for side in 1:2
    # Align interface in positive coordinate direction of primary element.
    # For orientation == 1, the secondary element needs to be indexed backwards
    # relative to the interface.
    if side == 1 || orientation == 0
      # Forward indexing
      i = :i_forward
    else
      # Backward indexing
      i = :i_backward
    end

    if faces[side] == 0
      # Index face in negative x-direction
      interfaces.node_indices[side, interface_id] = (:begin, i)
    elseif faces[side] == 1
      # Index face in positive x-direction
      interfaces.node_indices[side, interface_id] = (:end, i)
    elseif faces[side] == 2
      # Index face in negative y-direction
      interfaces.node_indices[side, interface_id] = (i, :begin)
    else # faces[side] == 3
      # Index face in positive y-direction
      interfaces.node_indices[side, interface_id] = (i, :end)
    end
  end

  return interfaces
end

# Initialize node_indices of boundary container.
@inline function init_boundary_node_indices!(boundaries::T8codeBoundaryContainer{2},
                                             face, boundary_id)
  if face == 0
    # Index face in negative x-direction
    boundaries.node_indices[boundary_id] = (:begin, :i_forward)
  elseif face == 1
    # Index face in positive x-direction
    boundaries.node_indices[boundary_id] = (:end, :i_forward)
  elseif face == 2
    # Index face in negative y-direction
    boundaries.node_indices[boundary_id] = (:i_forward, :begin)
  else # face == 3
    # Index face in positive y-direction
    boundaries.node_indices[boundary_id] = (:i_forward, :end)
  end

  return boundaries
end

# Initialize node_indices of mortar container
# faces[1] is expected to be the face of the small side.
@inline function init_mortar_node_indices!(mortars, faces, orientation, mortar_id)
  for side in 1:2
    # Align mortar in positive coordinate direction of small side.
    # For orientation == 1, the large side needs to be indexed backwards
    # relative to the mortar.
    if side == 1 || orientation == 0
      # Forward indexing for small side or orientation == 0
      i = :i_forward
    else
      # Backward indexing for large side with reversed orientation
      i = :i_backward
    end

    if faces[side] == 0
      # Index face in negative x-direction
      mortars.node_indices[side, mortar_id] = (:begin, i)
    elseif faces[side] == 1
      # Index face in positive x-direction
      mortars.node_indices[side, mortar_id] = (:end, i)
    elseif faces[side] == 2
      # Index face in negative y-direction
      mortars.node_indices[side, mortar_id] = (i, :begin)
    else # faces[side] == 3
      # Index face in positive y-direction
      mortars.node_indices[side, mortar_id] = (i, :end)
    end
  end

  return mortars
end

end # @muladd
