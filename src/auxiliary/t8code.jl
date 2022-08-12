T8DIR = ENV["JULIA_T8CODE_PATH"]

using Printf

# foobar
# println("aux = ", ENV["JULIA_T8CODE_PATH"])

libt8 = "$(T8DIR)/lib/libt8.so"
# libsc = "$(T8DIR)/lib/libsc.so"
# libp4 = "/home/jmark/install/t8code/lib/libt8.so"

struct t8_forest end

Cptr = Ptr{Cvoid}
CChar = UInt8
MPI_Comm_t = Cptr
t8_locidx_t = Int32

const SC_LP_DEFAULT    =  -1     # /**< this selects the SC default threshold */
const SC_LP_ALWAYS     =   0     # /**< this will log everything */
const SC_LP_TRACE      =   1     # /**< this will prefix file and line number */
const SC_LP_DEBUG      =   2     # /**< any information on the internal state */
const SC_LP_VERBOSE    =   3     # /**< information on conditions, decisions */
const SC_LP_INFO       =   4     # /**< the main things a function is doing */
const SC_LP_STATISTICS =   5     # /**< important for consistency/performance */
const SC_LP_PRODUCTION =   6     # /**< a few lines for a major api function */
const SC_LP_ESSENTIAL  =   7     # /**< this logs a few lines max per program */
const SC_LP_ERROR      =   8     # /**< this logs errors only */
const SC_LP_SILENT     =   9     # /**< this never logs anything */

@enum t8_vtk_data_type_t begin
  T8_VTK_SCALAR                # /* One double value per element */
  T8_VTK_VECTOR                # /* 3 double values per element */
end

macro SC_ASSERT(q)
  :( $(esc(q)) ? nothing : throw(AssertionError($(string(q)))) )
end

macro P4EST_ASSERT(q)
  :( @SC_ASSERT($(esc(q)) ) )
end

# P4EST_QUADRANT_LEN(l) = 1 << (P4EST_MAXLEVEL - l)

macro T8_ASSERT(q)
  :( @SC_ASSERT($(esc(q)) ) )
end

# This macro is not supposed to be understood by weaklings ...
macro t8_ccall(fname, rtype, args...)
  quote
    function $(esc(fname))($((typeof(arg.args[1]) == Symbol ? esc(arg.args[1]) : esc(Expr(:kw, arg.args[1].args[1], arg.args[2])) for arg in args)...))
      ccall(($(string(fname)), libt8), 
        $(esc(rtype)), ($((esc(typeof(arg.args[1]) == Symbol ? arg.args[2] : arg.args[1].args[2]) for arg in args)...),), 
        $((esc(typeof(arg.args[1]) == Symbol ? arg.args[1] : arg.args[1].args[1]) for arg in args)...))
    end
  end
end

# typedef int         (*t8_forest_adapt_t) (t8_forest_t forest,
#                                           t8_forest_t forest_from,
#                                           t8_locidx_t which_tree,
#                                           t8_locidx_t lelement_id,
#                                           t8_eclass_scheme_c *ts,
#                                           const int is_family,
#                                           const int num_elements,
#                                           t8_element_t *elements[]);
macro t8_adapt_callback(callback)
  :( @cfunction($callback, Cint, (Cptr, Cptr, t8_locidx_t, t8_locidx_t, Cptr, Cint, Cint, Ptr{Cptr})) )
end

@inline t8_mpi_comm() = mpi_comm().val

const T8CODE_MAXLEVEL = 30

# Macros from `t8code`
const t8code_root_len = 1 << T8CODE_MAXLEVEL
@inline t8code_quadrant_len(l) = 1 << (T8CODE_MAXLEVEL - l)

@t8_ccall(sc_init, Cvoid, comm :: MPI_Comm_t, catch_signals :: Cint, print_backtrace :: Cint, log_handler :: Cptr, log_threshold :: Cint)
@t8_ccall(t8_init, Cvoid, log_threshold :: Cint = SC_LP_PRODUCTION)

function init_t8code()
  # loglevel = SC_LP_VERBOSE
  # loglevel = SC_LP_SILENT
  loglevel = SC_LP_PRODUCTION

  sc_init(t8_mpi_comm(), 1, 1, C_NULL, loglevel)
  t8_init(loglevel)
  return nothing
end

function finalize_t8code()
  sc_finalize();
  return nothing
end

# int t8_get_package_id (void);
@t8_ccall(t8_get_package_id, Cint)

# void sc_free (int package, void *ptr);
@t8_ccall(sc_free, Cvoid, package :: Cint, ptr :: Ptr{Cvoid})

# int sc_memory_status (int package);
@t8_ccall(sc_memory_status, Cint, package :: Cint)

# void sc_memory_check (int package);
@t8_ccall(sc_memory_check, Cvoid, package :: Cint)

function t8_free(ptr)
  # sc_free(t8_get_package_id(), ptr)
  sc_free(-1, ptr)
end

@t8_ccall(t8_forest_init, Cvoid, pforest :: Cptr)

# void t8_forest_set_user_data (t8_forest_t forest, void *data);
@t8_ccall(t8_forest_set_user_data, Cvoid, pforest :: Cptr, data :: Cptr)

# void *t8_forest_get_user_data (t8_forest_t forest);
@t8_ccall(t8_forest_get_user_data, Cptr, pforest :: Cptr)

# int t8_forest_is_committed (t8_forest_t forest);
@t8_ccall(t8_forest_is_committed, Cint, pforest :: Cptr)

# void t8_forest_set_adapt (t8_forest_t forest,
#                           const t8_forest_t set_from,
#                           t8_forest_adapt_t adapt_fn,
#                           int recursive);
@t8_ccall(t8_forest_set_adapt, Cvoid, pforest :: Cptr, set_from :: Cptr, adapt_fn :: Cptr, recursive :: Cint)

# void t8_forest_set_partition (t8_forest_t forest,
#                               const t8_forest_t set_from,
#                               int set_for_coarsening);
@t8_ccall(t8_forest_set_partition, Cvoid, pforest :: Cptr, set_from :: Cptr, set_for_coarsening :: Cint)

# void t8_forest_set_balance (t8_forest_t forest,
#                             const t8_forest_t set_from,
#                             int no_repartition);
@t8_ccall(t8_forest_set_balance, Cvoid, pforest :: Cptr, set_from :: Cptr, no_repartition :: Cint)

# void t8_forest_commit (t8_forest_t forest);
@t8_ccall(t8_forest_commit, Cvoid, pforest :: Cptr)

# t8_locidx_t t8_forest_get_tree_element_offset (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_get_tree_element_offset, t8_locidx_t, forest :: Cptr, ltreeid :: t8_locidx_t)

@t8_ccall(t8_cmesh_new_periodic, Cptr, comm :: MPI_Comm_t, ndim :: Cint)
@t8_ccall(t8_cmesh_new_periodic_hybrid, Cptr, comm :: MPI_Comm_t)
@t8_ccall(t8_forest_get_num_global_trees, t8_locidx_t, forest :: Cptr)
@t8_ccall(t8_forest_get_global_num_elements, t8_locidx_t, forest :: Cptr)
@t8_ccall(t8_cmesh_new_from_p4est, Cptr, conn :: Cptr, comm :: MPI_Comm_t, do_partition :: Cint)

# t8_locidx_t t8_forest_get_num_ghosts (t8_forest_t forest);
@t8_ccall(t8_forest_get_num_ghosts, t8_locidx_t, forest :: Cptr)

@t8_ccall(t8_forest_get_num_local_trees, t8_locidx_t, forest :: Cptr)

@t8_ccall(t8_forest_get_local_num_elements, t8_locidx_t, forest :: Cptr)

# t8_locidx_t t8_forest_get_tree_num_elements (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_get_tree_num_elements , t8_locidx_t, forest :: Cptr, ltreeid :: t8_locidx_t)

# t8_eclass_t t8_forest_get_tree_class (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_get_tree_class, Cptr, forest :: Cptr, ltreeid :: t8_locidx_t)

# t8_eclass_scheme_c *t8_forest_get_eclass_scheme (t8_forest_t forest, t8_eclass_t eclass);
@t8_ccall(t8_forest_get_eclass_scheme, Cptr, forest :: Cptr, eclass :: Cptr)

# void t8_forest_element_centroid (t8_forest_t forest,
#                                  t8_locidx_t ltreeid,
#                                  const t8_element_t *element,
#                                  double *coordinates);
@t8_ccall(t8_forest_element_centroid, Cvoid, forest :: Cptr, ltreeid :: t8_locidx_t, element :: Cptr, coordinates :: Ptr{Cdouble})

# t8_element_t       *t8_forest_get_element_in_tree (t8_forest_t forest,
#                                                    t8_locidx_t ltreeid,
#                                                    t8_locidx_t leid_in_tree);
@t8_ccall(t8_forest_get_element_in_tree, Cptr, forest :: Cptr, ltreeid :: t8_locidx_t, leid_in_tree :: t8_locidx_t)



# function t8_cmesh_new_from_p4est(conn, do_partition = 0) :: Cptr
#   # cmesh = c"t8_cmesh_new_from_p4est"(conn, comm, do_partition)
#   # cmesh = C_NULL
#   # return cmesh
#   return C_NULL
# 
#   return ccall(("t8_cmesh_new_from_p4est", libt8), Cptr, (Cptr, MPI_Comm_t, Cint), conn, comm, do_partition)
# end

# function t8_scheme_new_default() :: Cptr
#   # t8_scheme_cxx_t    *t8_scheme_new_default_cxx (void);
#   return ccall(("t8_scheme_new_default_cxx", libt8), Cptr, ())
# end

@t8_ccall(t8_scheme_new_default_cxx, Cptr)

# function t8_forest_new_uniform(cmesh, scheme, level) :: Cptr
#   # t8_forest_t         t8_forest_new_uniform (t8_cmesh_t cmesh,
#   #                                         t8_scheme_cxx_t *scheme,
#   #                                         int level, int do_face_ghost,
#   #                                         sc_MPI_Comm comm);
#   return C_NULL
# 
#   do_face_ghost = 0
#   return ccall(("t8_forest_new_uniform", libt8), Cptr, 
#     (Cptr, Cptr, Cint, Cint, MPI_Comm_t), cmesh, scheme, level, do_face_ghost, comm)
# end

@t8_ccall(t8_forest_new_uniform, Cptr, cmesh :: Cptr, scheme :: Cptr, level :: Cint, do_face_ghost :: Cint, comm :: MPI_Comm_t)

# void t8_forest_unref (t8_forest_t *pforest);
@t8_ccall(t8_forest_unref, Cvoid, pforest :: Cptr)


# int t8_element_level (t8_eclass_scheme_c *ts, const t8_element_t *elem);
@t8_ccall(t8_element_level, Cint, ts :: Cptr, elem :: Cptr)

# double              t8_forest_element_volume (t8_forest_t forest,
#                                               t8_locidx_t ltreeid,
#                                               const t8_element_t *element);
@t8_ccall(t8_forest_element_volume, Cdouble, forest :: Cptr, ltreeid :: t8_locidx_t, element :: Cptr)

# /** Compute the coordinates of a given vertex of an element if a geometry
#  * for this tree is registered in the forest's cmesh.
#  * \param [in]      forest     The forest.
#  * \param [in]      ltree_id   The forest local id of the tree in which the element is.
#  * \param [in]      element    The element.
#  * \param [in]      corner_number The corner number, in Z-order, of the vertex which should be computed.
#  * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
#  *                             the x, y and z coordinates of the vertex.
#  */
# void                t8_forest_element_coordinate (t8_forest_t forest,
#                                                   t8_locidx_t ltree_id,
#                                                   const t8_element_t *element,
#                                                   int corner_number,
#                                                   double *coordinates);
@t8_ccall(t8_forest_element_coordinate, Cvoid, 
                                        forest        :: Cptr,
                                        ltreeid       :: t8_locidx_t,
                                        element       :: Cptr,
                                        corner_number :: Cint,
                                        coordinates   :: Ptr{Cdouble})

# void t8_forest_element_centroid (t8_forest_t forest,
#                                  t8_locidx_t ltreeid,
#                                  const t8_element_t *element,
#                                  double *coordinates);
@t8_ccall(t8_forest_element_centroid, Cvoid, forest :: Cptr, ltreeid :: t8_locidx_t, element :: Cptr, coordinates :: Ptr{Cdouble})

# /** Compute the number of corners of a given element.
#   * \param [in] elem The element.
#   * \return          The number of corners of \a elem.
#   */
# int                 t8_element_num_corners (const t8_element_t *elem);
@t8_ccall(t8_element_num_corners, Cint, ts :: Cptr, elem :: Cptr)

# /** Compute the number of faces of a given element.
#  * \param [in] elem The element.
#  * \return          The number of faces of \a elem.
#  */
# int                 t8_element_num_faces (const t8_element_t *elem);
@t8_ccall(t8_element_num_faces, Cint, ts :: Cptr, elem :: Cptr)

# /** Compute the leaf face neighbors of a forest.
#  * \param [in]    forest  The forest. Must have a valid ghost layer.
#  * \param [in]    ltreeid A local tree id.
#  * \param [in]    leaf    A leaf in tree \a ltreeid of \a forest.
#  * \param [out]   neighbor_leafs Unallocated on input. On output the neighbor
#  *                        leafs are stored here.
#  * \param [in]    face    The index of the face across which the face neighbors
#  *                        are searched.
#  * \param [out]   dual_face On output the face id's of the neighboring elements' faces.
#  * \param [out]   num_neighbors On output the number of neighbor leafs.
#  * \param [out]   pelement_indices Unallocated on input. On output the element indices
#  *                        of the neighbor leafs are stored here.
#  *                        0, 1, ... num_local_el - 1 for local leafs and
#  *                        num_local_el , ... , num_local_el + num_ghosts - 1 for ghosts.
#  * \param [out]   pneigh_scheme On output the eclass scheme of the neighbor elements.
#  * \param [in]    forest_is_balanced True if we know that \a forest is balanced, false
#  *                        otherwise.
#  * \note If there are no face neighbors, then *neighbor_leafs = NULL, num_neighbors = 0,
#  * and *pelement_indices = NULL on output.
#  * \note Currently \a forest must be balanced.
#  * \note \a forest must be committed before calling this function.
#  */
# void                t8_forest_leaf_face_neighbors (t8_forest_t forest,
#                                                    t8_locidx_t ltreeid,
#                                                    const t8_element_t *leaf,
#                                                    t8_element_t
#                                                    **pneighbor_leafs[],
#                                                    int face,
#                                                    int *dual_faces[],
#                                                    int *num_neighbors,
#                                                    t8_locidx_t
#                                                    **pelement_indices,
#                                                    t8_eclass_scheme_c
#                                                    **pneigh_scheme,
#                                                    int forest_is_balanced);
@t8_ccall(t8_forest_leaf_face_neighbors, Cvoid,
  forest              :: Cptr,
  ltreeid             :: t8_locidx_t,
  leaf                :: Cptr,
  pneighbor_leafs     :: Cptr,
  face                :: Cint,
  dual_faces          :: Cptr,
  num_neighbors       :: Cptr,
  pelement_indices    :: Cptr,
  pneigh_scheme       :: Cptr,
  forest_is_balanced  :: Cint)

function trixi_t8_unref_forest(forest :: Cptr)
  t8_forest_unref(Ref(forest))
end

function trixi_t8_count_interfaces(forest :: Cptr)
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(forest) != 0);

  # /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements(forest);
  # /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts(forest);
  # /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees(forest);

  current_index = 0

  local_num_conform = 0
  local_num_mortars = 0
  local_num_boundry = 0

  for itree = 0:num_local_trees-1
    tree_class = t8_forest_get_tree_class(forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    # /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      element = t8_forest_get_element_in_tree(forest, itree, ielement);

      level = t8_element_level(eclass_scheme, element)

      num_faces = t8_element_num_faces(eclass_scheme,element)

      for iface = 0:num_faces-1

        pelement_indices_ref = Ref{Ptr{t8_locidx_t}}()
        pneighbor_leafs_ref = Ref{Ptr{Cptr}}()
        pneigh_scheme_ref = Ref{Cptr}()

        dual_faces_ref = Ref{Ptr{Cint}}()
        num_neighbors_ref = Ref{Cint}()

        forest_is_balanced = Cint(1)

        t8_forest_leaf_face_neighbors(forest,itree,element,
          pneighbor_leafs_ref, iface, dual_faces_ref, num_neighbors_ref,
          pelement_indices_ref, pneigh_scheme_ref, forest_is_balanced)

        num_neighbors      = num_neighbors_ref[]
        dual_faces         = unsafe_wrap(Array,dual_faces_ref[],num_neighbors)
        neighbor_ielements = unsafe_wrap(Array,pelement_indices_ref[],num_neighbors)
        neighbor_leafs     = unsafe_wrap(Array,pneighbor_leafs_ref[],num_neighbors)
        neighbor_scheme    = pneigh_scheme_ref[]

        if num_neighbors > 0
          neighbor_level = t8_element_level(neighbor_scheme, neighbor_leafs[1])

          if level < neighbor_level 
              local_num_mortars += 1
          elseif level == neighbor_level && all(Int32(current_index) .<= neighbor_ielements)
              local_num_conform += 1
          end

        else

          local_num_boundry += 1

        end
       
        t8_free(dual_faces_ref[])
        t8_free(pneighbor_leafs_ref[])
        t8_free(pelement_indices_ref[])

      end # for

      current_index += 1
    end # for
  end # for

  # println("")
  # println("")
  # println(" ## local_num_elements = ", num_local_elements)
  # println(" ## local_num_conform  = ", local_num_conform)
  # println(" ## local_num_mortars  = ", local_num_mortars)
  # println(" ## local_num_boundry  = ", local_num_boundry)
  # println("")
  # println("")

  return (interfaces = local_num_conform,
          mortars    = local_num_mortars,
          boundaries = local_num_boundry)
end

function trixi_t8_fill_interfaces(forest :: Cptr, interfaces, mortars, boundaries)
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(forest) != 0);

  # /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements(forest);
  # /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts(forest);
  # /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees(forest);

  current_index = 0

  local_num_conform = 0
  local_num_mortars = 0
  local_num_boundry = 0

  for itree = 0:num_local_trees-1
    tree_class = t8_forest_get_tree_class(forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    # /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      element = t8_forest_get_element_in_tree(forest, itree, ielement);

      level = t8_element_level(eclass_scheme, element)

      num_faces = t8_element_num_faces(eclass_scheme,element)

      for iface = 0:num_faces-1

        pelement_indices_ref = Ref{Ptr{t8_locidx_t}}()
        pneighbor_leafs_ref = Ref{Ptr{Cptr}}()
        pneigh_scheme_ref = Ref{Cptr}()

        dual_faces_ref = Ref{Ptr{Cint}}()
        num_neighbors_ref = Ref{Cint}()

        forest_is_balanced = Cint(1)

        t8_forest_leaf_face_neighbors(forest,itree,element,
          pneighbor_leafs_ref, iface, dual_faces_ref, num_neighbors_ref,
          pelement_indices_ref, pneigh_scheme_ref, forest_is_balanced)

        num_neighbors      = num_neighbors_ref[]
        dual_faces         = unsafe_wrap(Array,dual_faces_ref[],num_neighbors)
        neighbor_ielements = unsafe_wrap(Array,pelement_indices_ref[],num_neighbors)
        neighbor_leafs     = unsafe_wrap(Array,pneighbor_leafs_ref[],num_neighbors)
        neighbor_scheme    = pneigh_scheme_ref[]

        if num_neighbors > 0
          neighbor_level = t8_element_level(neighbor_scheme, neighbor_leafs[1])

          # Non-conforming interface.
          if level < neighbor_level 
              local_num_mortars += 1

              # faces = (iface,dual_faces[1])
              faces = (dual_faces[1],iface)

              mortar_id = local_num_mortars
              orientation = 0 # aligned in z-order

              # Last entry is the large element ... What a stupid convention!
              mortars.neighbor_ids[end, mortar_id] = ielement + 1

              # First `1:end-1` entries are the smaller elements.
              mortars.neighbor_ids[1:end-1, mortar_id] .= neighbor_ielements[:] .+ 1

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

          # Conforming interface.
          elseif level == neighbor_level && all(Int32(current_index) .<= neighbor_ielements)
              local_num_conform += 1

              faces = (iface, dual_faces[1])
              interface_id = local_num_conform
              orientation = 0 # aligned in z-order

              # Write data to interfaces container.
              interfaces.neighbor_ids[1, interface_id] = ielement + 1
              interfaces.neighbor_ids[2, interface_id] = neighbor_ielements[1] + 1

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
          end

        # Domain boundary.
        else

          local_num_boundry += 1

        end
     
        t8_free(dual_faces_ref[])
        t8_free(pneighbor_leafs_ref[])
        t8_free(pelement_indices_ref[])

      end # for

      current_index += 1
    end # for
  end # for

  # println("")
  # println("")
  # println(" ## local_num_conform = ", local_num_conform)
  # println(" ## local_num_mortars = ", local_num_mortars)
  # println(" ## local_num_boundry = ", local_num_boundry)
  # println("")
  # println("")

  return (interfaces = local_num_conform,
          mortars    = local_num_mortars,
          boundaries = local_num_boundry)
end

function trixi_t8_get_local_element_levels(forest :: Cptr)
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(forest) != 0);

  levels = Vector{Int}(undef, t8_forest_get_local_num_elements(forest))

  # /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees(forest);

  current_index = 0

  for itree = 0:num_local_trees-1
    tree_class = t8_forest_get_tree_class(forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    # /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      element = t8_forest_get_element_in_tree(forest, itree, ielement);
      current_index += 1
      levels[current_index] = t8_element_level(eclass_scheme, element)
    end # for
  end # for

  return levels
end

function adapt_callback(forest,
                         forest_from,
                         which_tree,
                         lelement_id,
                         ts,
                         is_family, 
                         num_elements,
                         elements) :: Cint

  num_levels = t8_forest_get_local_num_elements(forest_from)

  indicator_ptr = Ptr{Int}(t8_forest_get_user_data(forest))
  indicators = unsafe_wrap(Array,indicator_ptr,num_levels)

  offset = t8_forest_get_tree_element_offset(forest_from, which_tree)

  # Only allow coarsening for complete families.
  if indicators[offset + lelement_id + 1] < 0 && is_family == 0
    return Cint(0)
  end

  return Cint(indicators[offset + lelement_id + 1])

end

function trixi_t8_adapt_new(old_forest :: Cptr, indicators)
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(old_forest) != 0);

  # Init new forest.
  new_forest_ref = Ref{Ptr{Cptr}}()
  t8_forest_init(new_forest_ref);
  new_forest = new_forest_ref[]

  let set_from = C_NULL, recursive = 0, set_for_coarsening = 0, no_repartition = 0
    t8_forest_set_user_data(new_forest, pointer(indicators))
    t8_forest_set_adapt(new_forest, old_forest, @t8_adapt_callback(adapt_callback), recursive)
    t8_forest_set_balance(new_forest, set_from, no_repartition)
    t8_forest_set_partition(new_forest, set_from, set_for_coarsening)
    # c"t8_forest_set_ghost(new_forest, 1, c"T8_GHOST_FACES");
    t8_forest_commit(new_forest)
  end

  return new_forest
end

function trixi_t8_get_difference(old_levels, new_levels)

  old_nelems = length(old_levels)
  new_nelems = length(new_levels)

  changes = Vector{Int}(undef, old_nelems)

  # Local element indices.
  old_index = 1
  new_index = 1

  # TODO: Make general for 2D/3D and hybrid grids.
  T8_CHILDREN = 4

  while old_index <= old_nelems && new_index <= new_nelems

    if old_levels[old_index] < new_levels[new_index] 
      # Refined.
    
      changes[old_index] = 1

      old_index += 1
      new_index += T8_CHILDREN

    elseif old_levels[old_index] > new_levels[new_index] 
      # Coarsend.

      for child_index = old_index:old_index+T8_CHILDREN-1
        changes[child_index] = -1
      end

      old_index += T8_CHILDREN
      new_index += 1

    else
      # No changes.
      
      changes[old_index] = 0

      old_index += 1
      new_index += 1
    end
  end

  return changes
end
