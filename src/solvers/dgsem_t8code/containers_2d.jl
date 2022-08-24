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
      element = t8_forest_get_element_in_tree(mesh.forest, itree, ielement)
      element_level  = t8_element_level(eclass_scheme, element)

      # TODO: Make this more general.
      element_length = t8code_quadrant_len(element_level) / t8code_root_len

      element_coords = Array{Float64}(undef,3)
      t8_element_vertex_reference_coords(eclass_scheme, element, 0, pointer(element_coords))

      nodes_out_x = 2 * (element_length * 1/2 * (nodes .+ 1) .+ element_coords[1]) .- 1
      nodes_out_y = 2 * (element_length * 1/2 * (nodes .+ 1) .+ element_coords[2]) .- 1

      polynomial_interpolation_matrix!(matrix1, mesh.nodes, nodes_out_x, baryweights_in)
      polynomial_interpolation_matrix!(matrix2, mesh.nodes, nodes_out_y, baryweights_in)

      multiply_dimensionwise!(
        view(node_coordinates, :, :, :, current_index += 1),
        matrix1, matrix2,
        view(mesh.tree_node_coordinates, :, :, :, itree+1),
        tmp1
      )
      # end TODO
    end
  end

  return node_coordinates
end

end # @muladd
