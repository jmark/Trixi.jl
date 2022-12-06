T8DIR = ENV["JULIA_T8CODE_PATH"]
libt8 = "$(T8DIR)/lib/libt8.so"

# libsc = "$(T8DIR)/lib/libsc.so"
# libp4 = "/home/jmark/install/t8code/lib/libt8.so"

const MPI_Comm_t = Ptr{Cvoid}
const t8_locidx_t = Int32
const t8_gloidx_t = Int64
const p4_topidx_t = Int32

const p4est_topidx_t = Int32
const p4est_gloidx_t = Int64

const p4est_t = Ptr{Cvoid}
const p8est_t = Ptr{Cvoid}

const p4est_quadrant_t = Ptr{Cvoid}
const p8est_quadrant_t = Ptr{Cvoid}

const p4est_iter_volume_info_t = Ptr{Cvoid}
const p8est_iter_volume_info_t = Ptr{Cvoid}

const p4est_iter_face_info_t = Ptr{Cvoid}
const p8est_iter_face_info_t = Ptr{Cvoid}

const p4est_iter_corner_info_t = Ptr{Cvoid}
const p8est_iter_corner_info_t = Ptr{Cvoid}

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

# /** The maximum number of boundary faces an element class can have. */
const T8_ECLASS_MAX_FACES = 6
# /** The maximum number of boundary edges an element class can have. */
const T8_ECLASS_MAX_EDGES = 12
# /** The maximum number of boundary edges a 2D element class can have. */
const T8_ECLASS_MAX_EDGES_2D = 4
# /** The maximum number of cornes a 2-dimensional element class can have. */
const T8_ECLASS_MAX_CORNERS_2D = 4
# /** The maximum number of cornes an element class can have. */
const T8_ECLASS_MAX_CORNERS = 8
# /** The maximal possible dimension for an eclass */
const T8_ECLASS_MAX_DIM = 3

@enum t8_vtk_data_type_t begin
  T8_VTK_SCALAR                # /* One double value per element */
  T8_VTK_VECTOR                # /* 3 double values per element */
end

# /** This enumeration contains all possible element classes. */
@enum t8_eclass_t begin
  # /** The vertex is the only zero-dimensional element class. */
  T8_ECLASS_VERTEX = 0
  # /** The line is the only one-dimensional element class. */
  T8_ECLASS_LINE
  # /** The quadrilateral is one of two element classes in two dimensions. */
  T8_ECLASS_QUAD
  # /** The element class for a triangle. */
  T8_ECLASS_TRIANGLE
  # /** The hexahedron is one three-dimensional element class. */
  T8_ECLASS_HEX
  # /** The tetrahedron is another three-dimensional element class. */
  T8_ECLASS_TET
  # /** The prism has five sides: two opposing triangles joined by three quadrilaterals. */
  T8_ECLASS_PRISM
  # /** The pyramid has a quadrilateral as base and four triangles as sides. */
  T8_ECLASS_PYRAMID
  # /** This is no element class but can be used as the number of element classes. */
  T8_ECLASS_COUNT
  # /** This is no element class but can be used for the case a class of a third party library is not supported by t8code*/
  T8_ECLASS_INVALID
end

t8_eclass_to_element_type = Dict(
  Int(T8_ECLASS_LINE)       => :line,
  Int(T8_ECLASS_QUAD)       => :quad,
  Int(T8_ECLASS_TRIANGLE)   => :tri,
  Int(T8_ECLASS_HEX)        => :hex,
  Int(T8_ECLASS_TET)        => :tet,
  Int(T8_ECLASS_PRISM)      => :prism,
  Int(T8_ECLASS_PYRAMID)    => :pyra,
)

t8_eclass_num_vertices = Dict(
  :line   => 2,
  :quad   => 4,
  :tri    => 3,
  :hex    => 8,
  :tet    => 4,
  :prism  => 6,
  :pyra   => 5,
)

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
  :( @cfunction($callback, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, t8_locidx_t, t8_locidx_t, Ptr{Cvoid}, Cint, Cint, Ptr{Ptr{Cvoid}})) )
end

@inline t8_mpi_comm() = mpi_comm().val

const T8CODE_MAXLEVEL = 30

# Macros from `t8code`
const t8code_root_len = 1 << T8CODE_MAXLEVEL
@inline t8code_quadrant_len(l) = 1 << (T8CODE_MAXLEVEL - l)

@t8_ccall(sc_init, Cvoid, comm :: MPI_Comm_t, catch_signals :: Cint, print_backtrace :: Cint, log_handler :: Ptr{Cvoid}, log_threshold :: Cint)

# /** Unregisters all packages, runs the memory check, removes the
#  * signal handlers and resets sc_identifier and sc_root_*.
#  * This function is optional.
#  * This function does not require sc_init to be called first.
#  */
# void sc_finalize (void);
@t8_ccall(sc_finalize, Cvoid)

@t8_ccall(t8_init, Cvoid, log_threshold :: Cint = SC_LP_PRODUCTION)

# void p4est_init (sc_log_handler_t log_handler, int log_threshold);
@t8_ccall(p4est_init, Cvoid, log_handler :: Ptr{Cvoid}, log_threshold :: Cint)

function init_t8code()
  # loglevel = SC_LP_VERBOSE
  loglevel = SC_LP_SILENT
  # loglevel = SC_LP_PRODUCTION

  sc_init(t8_mpi_comm(), 1, 1, C_NULL, loglevel)
  p4est_init(C_NULL, loglevel)
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

@t8_ccall(t8_forest_init, Cvoid, pforest :: Ptr{Cvoid})

# void t8_forest_set_user_data (t8_forest_t forest, void *data);
@t8_ccall(t8_forest_set_user_data, Cvoid, pforest :: Ptr{Cvoid}, data :: Ptr{Cvoid})

# void *t8_forest_get_user_data (t8_forest_t forest);
@t8_ccall(t8_forest_get_user_data, Ptr{Cvoid}, pforest :: Ptr{Cvoid})

# int t8_forest_is_committed (t8_forest_t forest);
@t8_ccall(t8_forest_is_committed, Cint, pforest :: Ptr{Cvoid})

# void t8_forest_set_adapt (t8_forest_t forest,
#                           const t8_forest_t set_from,
#                           t8_forest_adapt_t adapt_fn,
#                           int recursive);
@t8_ccall(t8_forest_set_adapt, Cvoid, pforest :: Ptr{Cvoid}, set_from :: Ptr{Cvoid}, adapt_fn :: Ptr{Cvoid}, recursive :: Cint)

# void t8_forest_set_partition (t8_forest_t forest,
#                               const t8_forest_t set_from,
#                               int set_for_coarsening);
@t8_ccall(t8_forest_set_partition, Cvoid, pforest :: Ptr{Cvoid}, set_from :: Ptr{Cvoid}, set_for_coarsening :: Cint)

# void t8_forest_set_balance (t8_forest_t forest,
#                             const t8_forest_t set_from,
#                             int no_repartition);
@t8_ccall(t8_forest_set_balance, Cvoid, pforest :: Ptr{Cvoid}, set_from :: Ptr{Cvoid}, no_repartition :: Cint)

# void t8_forest_commit (t8_forest_t forest);
@t8_ccall(t8_forest_commit, Cvoid, pforest :: Ptr{Cvoid})

# t8_locidx_t t8_forest_get_tree_element_offset (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_get_tree_element_offset, t8_locidx_t, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t)

@t8_ccall(t8_cmesh_new_periodic, Ptr{Cvoid}, comm :: MPI_Comm_t, ndim :: Cint)

@t8_ccall(t8_cmesh_new_periodic_hybrid, Ptr{Cvoid}, comm :: MPI_Comm_t)

@t8_ccall(t8_forest_get_num_global_trees, t8_locidx_t, forest :: Ptr{Cvoid})

@t8_ccall(t8_forest_get_global_num_elements, t8_locidx_t, forest :: Ptr{Cvoid})

# /** Return a cmesh associated to a forest.
#  * \param [in]      forest      The forest.
#  * \return          The cmesh associated to the forest.
#  */
# t8_cmesh_t          t8_forest_get_cmesh (t8_forest_t forest);
@t8_ccall(t8_forest_get_cmesh, Ptr{Cvoid}, forest :: Ptr{Cvoid})

# t8_cmesh_t t8_cmesh_new_from_p4est (p4est_connectivity_t * conn, sc_MPI_Comm comm, int do_partition);
@t8_ccall(t8_cmesh_new_from_p4est, Ptr{Cvoid}, conn :: Ptr{Cvoid}, comm :: MPI_Comm_t, do_partition :: Cint)

# t8_cmesh_t t8_cmesh_new_from_p8est (p8est_connectivity_t * conn, sc_MPI_Comm comm, int do_partition);
@t8_ccall(t8_cmesh_new_from_p8est, Ptr{Cvoid}, conn :: Ptr{Cvoid}, comm :: MPI_Comm_t, do_partition :: Cint)

# t8_cmesh_t t8_cmesh_from_msh_file (const char *fileprefix, int partition, sc_MPI_Comm comm, int dim, int master, int use_occ_geometry);
@t8_ccall(t8_cmesh_from_msh_file, Ptr{Cvoid}, fileprefix :: Cstring, partition :: Cint, comm :: MPI_Comm_t, dim :: Cint, master :: Cint, use_occ_geometry :: Cint)

# double *t8_cmesh_get_tree_vertices (t8_cmesh_t cmesh, t8_locidx_t ltreeid);
@t8_ccall(t8_cmesh_get_tree_vertices, Ptr{Cdouble}, cmesh :: Ptr{Cvoid}, ltreeid :: t8_locidx_t)

# /** Return the eclass of a given local tree.
#  * TODO: Should we refer to indices or consequently use ctree_t?
#  * \param [in]    cmesh         The cmesh to be considered.
#  * \param [in]    tree_id       The local id of the tree whose eclass will be returned.
#  * \return                      The eclass of the given tree.
#  * TODO: Call tree ids ltree_id or gtree_id etc. instead of tree_id.
#  * \a cmesh must be committed before calling this function.
#  */
# t8_eclass_t t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_locidx_t ltree_id);
@t8_ccall(t8_cmesh_get_tree_class, Cint, cmesh :: Ptr{Cvoid}, ltreeid :: t8_locidx_t)

# /** Return the global number of trees in a cmesh.
#  * \param [in] cmesh       The cmesh to be considered.
#  * \return                 The number of trees associated to \a cmesh.
#  * \a cmesh must be committed before calling this function.
#  */
# t8_gloidx_t         t8_cmesh_get_num_trees (t8_cmesh_t cmesh);
@t8_ccall(t8_cmesh_get_num_trees, t8_gloidx_t, cmesh :: Ptr{Cvoid})

# /** Return the number of local trees of a cmesh.
#  *  If the cmesh is not partitioned this is equivalent to \ref t8_cmesh_get_num_trees.
#  * \param [in] cmesh       The cmesh to be considered.
#  * \return                 The number of local trees of the cmesh.
#  * \a cmesh must be committed before calling this function.
#  */
# t8_locidx_t         t8_cmesh_get_num_local_trees (t8_cmesh_t cmesh);
@t8_ccall(t8_cmesh_get_num_local_trees, t8_locidx_t, cmesh :: Ptr{Cvoid})

# /** Given a local tree id and a face number, get information about the face neighbor tree.
#  * \param [in]      cmesh     The cmesh to be considered.
#  * \param [in]      ltreeid   The local id of a tree or a ghost.
#  * \param [in]      face      A face number of the tree/ghost.
#  * \param [out]     dual_face If not NULL, the face number of the neighbor tree at this connection.
#  * \param [out]     orientation If not NULL, the face orientation of the connection.
#  * \return                    If non-negative: The local id of the neighbor tree or ghost.
#  *                            If negative: There is no neighbor across this face. \a dual_face and
#  *                            \a orientation remain unchanged.
#  * \note If \a ltreeid is a ghost and it has a neighbor which is neither a local tree or ghost,
#  *       then the return value will be negative.
#  *       Thus, a negative return value does not necessarily mean that this is a domain boundary.
#  *       To find out whether a tree is a domain boundary or not \see t8_cmesh_tree_face_is_boundary.
#  */
# t8_locidx_t         t8_cmesh_get_face_neighbor (const t8_cmesh_t cmesh,
#                                                 const t8_locidx_t ltreeid,
#                                                 const int face,
#                                                 int *dual_face,
#                                                 int *orientation);
@t8_ccall(t8_cmesh_get_face_neighbor, t8_locidx_t, cmesh :: Ptr{Cvoid},
  ltreeid :: t8_locidx_t, face :: Cint, dual_face :: Ptr{Cint}, orientation :: Ptr{Cint})

# t8_locidx_t t8_forest_get_num_ghosts (t8_forest_t forest);
@t8_ccall(t8_forest_get_num_ghosts, t8_locidx_t, forest :: Ptr{Cvoid})

@t8_ccall(t8_forest_get_num_local_trees, t8_locidx_t, forest :: Ptr{Cvoid})

@t8_ccall(t8_forest_get_local_num_elements, t8_locidx_t, forest :: Ptr{Cvoid})

# t8_locidx_t t8_forest_get_tree_num_elements (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_get_tree_num_elements , t8_locidx_t, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t)

# t8_eclass_t t8_forest_get_tree_class (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_get_tree_class, Ptr{Cvoid}, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t)

# t8_eclass_scheme_c *t8_forest_get_eclass_scheme (t8_forest_t forest, t8_eclass_t eclass);
@t8_ccall(t8_forest_get_eclass_scheme, Ptr{Cvoid}, forest :: Ptr{Cvoid}, eclass :: Ptr{Cvoid})

# t8_element_t       *t8_forest_get_element_in_tree (t8_forest_t forest,
#                                                    t8_locidx_t ltreeid,
#                                                    t8_locidx_t leid_in_tree);
@t8_ccall(t8_forest_get_element_in_tree, Ptr{Cvoid}, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t, leid_in_tree :: t8_locidx_t)

# /** Compute the coordinates of a given element vertex inside a reference tree
#  *  that is embedded into [0,1]^d (d = dimension).
#  * \param [in] ts     Implementation of a class scheme.
#  * \param [in] t      The element to be considered.
#  * \param [in] vertex The id of the vertex whose coordinates shall be computed.
#  * \param [out] coords An array of at least as many doubles as the element's dimension
#  *                      whose entries will be filled with the coordinates of \a vertex.
#  */
# void t8_element_vertex_reference_coords (t8_eclass_scheme_c *ts, const t8_element_t *t, int vertex, double coords[]);
@t8_ccall(t8_element_vertex_reference_coords, Cvoid, ts :: Ptr{Cvoid}, element :: Ptr{Cvoid}, vertex :: Cint, coordinates :: Ptr{Cdouble})

# function t8_cmesh_new_from_p4est(conn, do_partition = 0) :: Ptr{Cvoid}
#   # cmesh = c"t8_cmesh_new_from_p4est"(conn, comm, do_partition)
#   # cmesh = C_NULL
#   # return cmesh
#   return C_NULL
# 
#   return ccall(("t8_cmesh_new_from_p4est", libt8), Ptr{Cvoid}, (Ptr{Cvoid}, MPI_Comm_t, Cint), conn, comm, do_partition)
# end

# function t8_scheme_new_default() :: Ptr{Cvoid}
#   # t8_scheme_cxx_t    *t8_scheme_new_default_cxx (void);
#   return ccall(("t8_scheme_new_default_cxx", libt8), Ptr{Cvoid}, ())
# end

@t8_ccall(t8_scheme_new_default_cxx, Ptr{Cvoid})

# /** Compute whether a given element shares a given face with its root tree.
#  * \param [in] ts       Implementation of a class scheme.
#  * \param [in] elem     The input element.
#  * \param [in] face     A face of \a elem.
#  * \return              True if \a face is a subface of the element's root element.
#  */
# int t8_element_is_root_boundary (t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);
@t8_ccall(t8_element_is_root_boundary, Cint, ts :: Ptr{Cvoid}, eleme :: Ptr{Cvoid}, face :: Cint)

# function t8_forest_new_uniform(cmesh, scheme, level) :: Ptr{Cvoid}
#   # t8_forest_t         t8_forest_new_uniform (t8_cmesh_t cmesh,
#   #                                         t8_scheme_cxx_t *scheme,
#   #                                         int level, int do_face_ghost,
#   #                                         sc_MPI_Comm comm);
#   return C_NULL
# 
#   do_face_ghost = 0
#   return ccall(("t8_forest_new_uniform", libt8), Ptr{Cvoid}, 
#     (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint, MPI_Comm_t), cmesh, scheme, level, do_face_ghost, comm)
# end

@t8_ccall(t8_forest_new_uniform, Ptr{Cvoid}, cmesh :: Ptr{Cvoid}, scheme :: Ptr{Cvoid}, level :: Cint, do_face_ghost :: Cint, comm :: MPI_Comm_t)

# void t8_forest_unref (t8_forest_t *pforest);
@t8_ccall(t8_forest_unref, Cvoid, pforest :: Ptr{Cvoid})


# int t8_element_level (t8_eclass_scheme_c *ts, const t8_element_t *elem);
@t8_ccall(t8_element_level, Cint, ts :: Ptr{Cvoid}, elem :: Ptr{Cvoid})

# /** Given an element and a face of this element. If the face lies on the
#  *  tree boundary, return the face number of the tree face.
#  *  If not the return value is arbitrary.
#  * \param [in] ts       Implementation of a class scheme.
#  * \param [in] elem     The element.
#  * \param [in] face     The index of a face of \a elem.
#  * \return The index of the tree face that \a face is a subface of, if
#  *         \a face is on a tree boundary.
#  *         Any arbitrary integer if \a is not at a tree boundary.
#  */
# int t8_element_tree_face (t8_eclass_scheme_c *ts, const t8_element_t *elem, int face);
@t8_ccall(t8_element_tree_face, Cint, ts :: Ptr{Cvoid}, elem :: Ptr{Cvoid}, face :: Cint)

# /** Return the shape of an allocated element according its type.
# *  For example, a child of an element can be an element of a different shape
# *  and has to be handled differently - according to its shape.
# *  \param [in] ts     Implementation of a class scheme.
# *  \param [in] elem   The element to be considered
# *  \return            The shape of the element as an eclass
# */
# t8_element_shape_t  t8_element_shape (t8_eclass_scheme_c *ts,
#                                       const t8_element_t *elem);
@t8_ccall(t8_element_shape, Cint, ts :: Ptr{Cvoid}, elem :: Ptr{Cvoid})

# /** Compute the shape of the face of an element.
#  * \param [in] ts       Implementation of a class scheme.
#  * \param [in] elem     The element.
#  * \param [in] face     A face of \a elem.
#  * \return              The element shape of the face.
#  * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
#  *      and depending on the face number either T8_ECLASS_QUAD or
#  *      T8_ECLASS_TRIANGLE for prisms.
#  */
# t8_element_shape_t t8_element_face_shape (t8_eclass_scheme_c *ts,
#                                            const t8_element_t *elem,
#                                            int face);
@t8_ccall(t8_element_face_shape, Cint, ts :: Ptr{Cvoid}, elem :: Ptr{Cvoid}, face :: Cint)

# double              t8_forest_element_volume (t8_forest_t forest,
#                                               t8_locidx_t ltreeid,
#                                               const t8_element_t *element);
@t8_ccall(t8_forest_element_volume, Cdouble, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t, element :: Ptr{Cvoid})

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
                                        forest        :: Ptr{Cvoid},
                                        ltreeid       :: t8_locidx_t,
                                        element       :: Ptr{Cvoid},
                                        corner_number :: Cint,
                                        coordinates   :: Ptr{Cdouble})

# /** Given the local id of a tree in a forest, compute the tree's local id
#  * in the associated cmesh.
#  *  \param [in] forest    The forest.
#  *  \param [in] ltreeid   The local id of a tree or ghost in the forest.
#  * \return  The local id of the tree in the cmesh associated with the forest.
#  * \a forest must be committed before calling this function.
#  * \note For forest local trees, this is the inverse function of \ref t8_forest_cmesh_ltreeid_to_ltreeid.
#  */
# t8_locidx_t t8_forest_ltreeid_to_cmesh_ltreeid (t8_forest_t forest, t8_locidx_t ltreeid);
@t8_ccall(t8_forest_ltreeid_to_cmesh_ltreeid, t8_locidx_t, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t)

# void t8_forest_element_centroid (t8_forest_t forest,
#                                  t8_locidx_t ltreeid,
#                                  const t8_element_t *element,
#                                  double *coordinates);
@t8_ccall(t8_forest_element_centroid, Cvoid, forest :: Ptr{Cvoid}, ltreeid :: t8_locidx_t, element :: Ptr{Cvoid}, coordinates :: Ptr{Cdouble})

# /** Compute the number of corners of a given element.
#   * \param [in] elem The element.
#   * \return          The number of corners of \a elem.
#   */
# int t8_element_num_corners (const t8_element_t *elem);
@t8_ccall(t8_element_num_corners, Cint, ts :: Ptr{Cvoid}, elem :: Ptr{Cvoid})

# /** Compute the number of faces of a given element.
#  * \param [in] elem The element.
#  * \return          The number of faces of \a elem.
#  */
# int t8_element_num_faces (const t8_element_t *elem);
@t8_ccall(t8_element_num_faces, Cint, ts :: Ptr{Cvoid}, elem :: Ptr{Cvoid})

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
  forest              :: Ptr{Cvoid},
  ltreeid             :: t8_locidx_t,
  leaf                :: Ptr{Cvoid},
  pneighbor_leafs     :: Ptr{Cvoid},
  face                :: Cint,
  dual_faces          :: Ptr{Cvoid},
  num_neighbors       :: Ptr{Cvoid},
  pelement_indices    :: Ptr{Cvoid},
  pneigh_scheme       :: Ptr{Cvoid},
  forest_is_balanced  :: Cint)

# /** Allocate a connectivity structure and populate from constants.
#  * The attribute fields are initialized to NULL.
#  * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
#  * \param [in] num_trees      Number of trees in the forest.
#  * \param [in] num_corners    Number of tree-connecting corners.
#  * \param [in] coff           Corner-to-tree offsets (num_corners + 1 values).
#  *                            This must always be non-NULL; in trivial cases
#  *                            it is just a pointer to a p4est_topix value of 0.
#  * \return                    The connectivity is checked for validity.
#  */
# p4est_connectivity_t *p4est_connectivity_new_copy (p4est_topidx_t num_vertices,
#                                                    p4est_topidx_t num_trees,
#                                                    p4est_topidx_t num_corners,
#                                                    const double *vertices,
#                                                    const p4est_topidx_t * ttv,
#                                                    const p4est_topidx_t * ttt,
#                                                    const int8_t * ttf,
#                                                    const p4est_topidx_t * ttc,
#                                                    const p4est_topidx_t * coff,
#                                                    const p4est_topidx_t * ctt,
#                                                    const int8_t * ctc);
@t8_ccall(p4est_connectivity_new_copy, Ptr{Cvoid},
  num_vertices  :: p4est_topidx_t,
  num_trees     :: p4est_topidx_t,
  num_corners   :: p4est_topidx_t,
  vertices      :: Ptr{Cdouble},
  ttv           :: Ptr{p4est_topidx_t},
  ttt           :: Ptr{p4est_topidx_t},
  ttf           :: Ptr{Int8},
  ttc           :: Ptr{p4est_topidx_t},
  coff          :: Ptr{p4est_topidx_t},
  ctt           :: Ptr{p4est_topidx_t},
  ctc           :: Ptr{Int8})

# /** Examine a connectivity structure.
#  * \return          Returns true if structure is valid, false otherwise.
#  */
# int p4est_connectivity_is_valid (p4est_connectivity_t *connectivity);
@t8_ccall(p4est_connectivity_is_valid, Cint, connectivity :: Ptr{Cvoid})

# void p4est_connectivity_destroy (p4est_connectivity_t *connectivity);
@t8_ccall(p4est_connectivity_destroy, Cvoid, connectivity :: Ptr{Cvoid})

# p4est_connectivity_t *p4est_connectivity_new_brick (int mi, int ni, int periodic_a, int periodic_b);
@t8_ccall(p4est_connectivity_new_brick, Ptr{Cvoid}, mi :: Cint, ni :: Cint, periodic_a :: Cint, periodic_b :: Cint)

# p4est_connectivity_t *p4est_connectivity_read_inp (const char *filename);
@t8_ccall(p4est_connectivity_read_inp, Ptr{Cvoid}, filename :: Cstring)

function trixi_t8_unref_forest(forest :: Ptr{Cvoid})
  t8_forest_unref(Ref(forest))
end

function trixi_t8_count_interfaces(forest :: Ptr{Cvoid})
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
        pneighbor_leafs_ref = Ref{Ptr{Ptr{Cvoid}}}()
        pneigh_scheme_ref = Ref{Ptr{Cvoid}}()

        dual_faces_ref = Ref{Ptr{Cint}}()
        num_neighbors_ref = Ref{Cint}()

        forest_is_balanced = Cint(1)

        t8_forest_leaf_face_neighbors(forest,itree,element,
          pneighbor_leafs_ref, iface, dual_faces_ref, num_neighbors_ref,
          pelement_indices_ref, pneigh_scheme_ref, forest_is_balanced)

        num_neighbors      = num_neighbors_ref[]
        neighbor_ielements = unsafe_wrap(Array,pelement_indices_ref[],num_neighbors)
        neighbor_leafs     = unsafe_wrap(Array,pneighbor_leafs_ref[],num_neighbors)
        neighbor_scheme    = pneigh_scheme_ref[]

        if num_neighbors > 0
          neighbor_level = t8_element_level(neighbor_scheme, neighbor_leafs[1])

          if level == neighbor_level && all(Int32(current_index) .<= neighbor_ielements)
          # TODO: Find a fix for the case: Single element on root level with periodic boundaries.
          # elseif level == neighbor_level && 
          #   (all(Int32(current_index) .< neighbor_ielements) || 
          #   level == 0 && (iface == 0 || iface == 2 || iface == 4))
              local_num_conform += 1
          elseif level < neighbor_level 
              local_num_mortars += 1
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
  # println(" ## local_num_elements = ", num_local_elements)
  # println(" ## local_num_conform  = ", local_num_conform)
  # println(" ## local_num_mortars  = ", local_num_mortars)
  # println(" ## local_num_boundry  = ", local_num_boundry)
  # println("")

  return (interfaces = local_num_conform,
          mortars    = local_num_mortars,
          boundaries = local_num_boundry)
end

function trixi_t8_fill_mesh_info(forest :: Ptr{Cvoid}, elements, interfaces, mortars, boundaries, boundary_names)
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

    # Get the number of elements of this tree.
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      element = t8_forest_get_element_in_tree(forest, itree, ielement);

      level = t8_element_level(eclass_scheme, element)

      element_shape = t8_eclass_to_element_type[t8_element_shape(eclass_scheme, element)]
      elements.shapes[current_index+1] = element_shape

      num_faces = t8_element_num_faces(eclass_scheme,element)

      for iface = 0:num_faces-1

        # Compute the `orientation` of the touching faces.
        if t8_element_is_root_boundary(eclass_scheme, element, iface) == 1
          cmesh = t8_forest_get_cmesh(forest)
          itree_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid(forest, itree)
          iface_in_tree = t8_element_tree_face(eclass_scheme, element, iface)
          orientation_ref = Ref{Cint}()

          t8_cmesh_get_face_neighbor(cmesh, itree_in_cmesh, iface_in_tree, C_NULL, orientation_ref)
          orientation = orientation_ref[]
        else
          orientation = 0
        end

        element_face_shape = t8_eclass_to_element_type[t8_element_face_shape(eclass_scheme, element, iface)]

        pelement_indices_ref = Ref{Ptr{t8_locidx_t}}()
        pneighbor_leafs_ref = Ref{Ptr{Ptr{Cvoid}}}()
        pneigh_scheme_ref = Ref{Ptr{Cvoid}}()

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

          # Conforming interface: The second condition ensures we only visit the interface once.
          if level == neighbor_level && Int32(current_index) <= neighbor_ielements[1]
          # elseif level == neighbor_level &&
          #   (all(Int32(current_index) .< neighbor_ielements) ||
          #   level == 0 && (iface == 0 || iface == 2 || iface == 4))
              local_num_conform += 1

              faces = (iface, dual_faces[1])
              interface_id = local_num_conform

              # Write data to interfaces container.
              interfaces.neighbor_ids[1, interface_id] = current_index + 1
              interfaces.neighbor_ids[2, interface_id] = neighbor_ielements[1] + 1

              # TODO: Is this correct for 3D, too?
              interfaces.shapes[interface_id] = element_face_shape

              # println("neighbor_ids = ", interfaces.neighbor_ids[:,interface_id] .- 1)
              # println("faces = ", faces)

              # Iterate over primary and secondary element.
              for side in 1:2
                # Align interface in positive coordinate direction of primary element.
                # For orientation == 1, the secondary element needs to be indexed backwards
                # relative to the interface.
                if side == 1 || orientation == 0
                  # Forward indexing
                  indexing = :i_forward
                else
                  # Backward indexing
                  indexing = :i_backward
                end

                if faces[side] == 0
                  # Index face in negative x-direction
                  interfaces.node_indices[side, interface_id] = (:begin, indexing)
                elseif faces[side] == 1
                  # Index face in positive x-direction
                  interfaces.node_indices[side, interface_id] = (:end, indexing)
                elseif faces[side] == 2
                  # Index face in negative y-direction
                  interfaces.node_indices[side, interface_id] = (indexing, :begin)
                else # faces[side] == 3
                  # Index face in positive y-direction
                  interfaces.node_indices[side, interface_id] = (indexing, :end)
                end
              end

          # Non-conforming interface.
          elseif level < neighbor_level 
              local_num_mortars += 1

              faces = (dual_faces[1],iface)

              mortar_id = local_num_mortars

              # Last entry is the large element ... What a stupid convention!
              # mortars.neighbor_ids[end, mortar_id] = ielement + 1
              mortars.neighbor_ids[end, mortar_id] = current_index + 1

              # First `1:end-1` entries are the smaller elements.
              mortars.neighbor_ids[1:end-1, mortar_id] .= neighbor_ielements[:] .+ 1

              # TODO: Is this correct for 3D, too?
              mortars.shapes[:, mortar_id] .= element_face_shape

              for side in 1:2
                # Align mortar in positive coordinate direction of small side.
                # For orientation == 1, the large side needs to be indexed backwards
                # relative to the mortar.
                if side == 1 || orientation == 0
                  # Forward indexing for small side or orientation == 0
                  indexing = :i_forward
                else
                  # Backward indexing for large side with reversed orientation
                  indexing = :i_backward
                end

                if faces[side] == 0
                  # Index face in negative x-direction
                  mortars.node_indices[side, mortar_id] = (:begin, indexing)
                elseif faces[side] == 1
                  # Index face in positive x-direction
                  mortars.node_indices[side, mortar_id] = (:end, indexing)
                elseif faces[side] == 2
                  # Index face in negative y-direction
                  mortars.node_indices[side, mortar_id] = (indexing, :begin)
                else # faces[side] == 3
                  # Index face in positive y-direction
                  mortars.node_indices[side, mortar_id] = (indexing, :end)
                end
              end
            
          # else: "level > neighbor_level" is skipped since we visit the mortar interface only once.
          end

        # Domain boundary.
        else
          local_num_boundry += 1
          boundary_id = local_num_boundry

          boundaries.neighbor_ids[boundary_id] = current_index + 1

          # println("neighbor_id = ", boundaries.neighbor_ids[boundary_id] -1)
          # println("face = ", iface)

          if iface == 0
            # Index face in negative x-direction.
            boundaries.node_indices[boundary_id] = (:begin, :i_forward)
          elseif iface == 1
            # Index face in positive x-direction.
            boundaries.node_indices[boundary_id] = (:end, :i_forward)
          elseif iface == 2
            # Index face in negative y-direction.
            boundaries.node_indices[boundary_id] = (:i_forward, :begin)
          else # iface == 3
            # Index face in positive y-direction.
            boundaries.node_indices[boundary_id] = (:i_forward, :end)
          end

          # One-based indexing.
          boundaries.name[boundary_id] = boundary_names[iface + 1, itree + 1]
        end
     
        t8_free(dual_faces_ref[])
        t8_free(pneighbor_leafs_ref[])
        t8_free(pelement_indices_ref[])

      end # for iface = ...

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

function trixi_t8_count_faces(forest :: Ptr{Cvoid})
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(forest) != 0);

  # /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements(forest);
  # /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts(forest);
  # /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees(forest);

  local_num_faces = 0

  for itree = 0:num_local_trees-1
    tree_class = t8_forest_get_tree_class(forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    # /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      element = t8_forest_get_element_in_tree(forest, itree, ielement)

      # local_num_faces += t8_element_num_faces(eclass_scheme,element)

      num_faces = t8_element_num_faces(eclass_scheme,element)

      for iface = 0:num_faces-1

        pelement_indices_ref = Ref{Ptr{t8_locidx_t}}()
        pneighbor_leafs_ref = Ref{Ptr{Ptr{Cvoid}}}()
        pneigh_scheme_ref = Ref{Ptr{Cvoid}}()

        dual_faces_ref = Ref{Ptr{Cint}}()
        num_neighbors_ref = Ref{Cint}()

        forest_is_balanced = Cint(1)

        t8_forest_leaf_face_neighbors(forest,itree,element,
          pneighbor_leafs_ref, iface, dual_faces_ref, num_neighbors_ref,
          pelement_indices_ref, pneigh_scheme_ref, forest_is_balanced)

        local_num_faces += max(1,num_neighbors_ref[])

        # println(max(1,num_neighbors_ref[]))
        # println(local_num_faces)

        # neighbor_ielements = unsafe_wrap(Array,pelement_indices_ref[],num_neighbors)
        # neighbor_leafs     = unsafe_wrap(Array,pneighbor_leafs_ref[],num_neighbors)
        # neighbor_scheme    = pneigh_scheme_ref[]

        # if num_neighbors > 0
        #   neighbor_level = t8_element_level(neighbor_scheme, neighbor_leafs[1])

        #   if level == neighbor_level && all(Int32(current_index) .<= neighbor_ielements)
        #   # TODO: Find a fix for the case: Single element on root level with periodic boundaries.
        #   # elseif level == neighbor_level && 
        #   #   (all(Int32(current_index) .< neighbor_ielements) || 
        #   #   level == 0 && (iface == 0 || iface == 2 || iface == 4))
        #       local_num_conform += 1
        #   elseif level < neighbor_level 
        #       local_num_mortars += 1
        #   end

        # else

        #   local_num_boundry += 1

        # end
       
        t8_free(dual_faces_ref[])
        t8_free(pneighbor_leafs_ref[])
        t8_free(pelement_indices_ref[])

      end # for

      # println(local_num_faces, " ", t8_element_num_faces(eclass_scheme,element))
    end # for
  end # for

  return local_num_faces

end

function trixi_t8_fill_mesh_info!(forest :: Ptr{Cvoid}, shapes, levels, mapP, mapM, orientation)
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(forest) != 0);

  # /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements(forest);
  # /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts(forest);
  # /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees(forest);

  elem_index = 0
  face_index = 0

  for itree = 0:num_local_trees-1
    tree_class = t8_forest_get_tree_class(forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    # Get the number of elements of this tree.
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1

      elem_index += 1

      element = t8_forest_get_element_in_tree(forest, itree, ielement);

      level = t8_element_level(eclass_scheme, element)

      levels[elem_index] = level
      shapes[elem_index] = t8_eclass_to_element_type[t8_element_shape(eclass_scheme, element)]

      num_faces = t8_element_num_faces(eclass_scheme,element)

      for iface = 0:num_faces-1

        # Compute the `orientation` of the touching faces.
        if t8_element_is_root_boundary(eclass_scheme, element, iface) == 1
          cmesh = t8_forest_get_cmesh(forest)
          itree_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid(forest, itree)
          iface_in_tree = t8_element_tree_face(eclass_scheme, element, iface)
          orientation_ref = Ref{Cint}()

          t8_cmesh_get_face_neighbor(cmesh, itree_in_cmesh, iface_in_tree, C_NULL, orientation_ref)
          orient = orientation_ref[]
        else
          orient = 0
        end

        # face_shapes[face_index] = t8_eclass_to_element_type[t8_element_face_shape(eclass_scheme, element, iface)]

        pelement_indices_ref = Ref{Ptr{t8_locidx_t}}()
        pneighbor_leafs_ref = Ref{Ptr{Ptr{Cvoid}}}()
        pneigh_scheme_ref = Ref{Ptr{Cvoid}}()

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

        # In case there is a boundary: `num_neighbors == 0`
        if num_neighbors < 1
          face_index += 1
          mapP[face_index] = elem_index
          mapM[face_index] = num_neighbors
          orientation[face_index] = orient
        else
          for ineighbor = 1:num_neighbors
            face_index += 1
            mapP[face_index + ineighbor-1] = elem_index
            mapM[face_index + ineighbor-1] = neighbor_ielements[ineighbor]
            orientation[face_index] = orient
          end
        end

        t8_free(dual_faces_ref[])
        t8_free(pneighbor_leafs_ref[])
        t8_free(pelement_indices_ref[])

      end # for iface = ...
    end # for
  end # for

end

function trixi_t8_get_local_element_levels(forest :: Ptr{Cvoid})
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

function trixi_t8_adapt_new(old_forest :: Ptr{Cvoid}, indicators)
  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(old_forest) != 0);

  # Init new forest.
  new_forest_ref = Ref{Ptr{Ptr{Cvoid}}}()
  t8_forest_init(new_forest_ref);
  new_forest :: Ptr{Cvoid} = new_forest_ref[]

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
