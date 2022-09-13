@muladd begin

# FIXME: This routine is identical to `dgsem_p4est/dg_2d.jl` but somehow there is an
# ambiguity between method definitions. This hinders to just use the implemenation from `dgsem_p4est`.
# The methods below are specialized on the mortar type
# and called from the basic `create_cache` method at the top.
function create_cache(mesh::T8codeMesh{2}, equations, mortar_l2::LobattoLegendreMortarL2, uEltype)
  # TODO: Total performance using different types.
  MA2d = MArray{Tuple{nvariables(equations), nnodes(mortar_l2)},
                uEltype, 2,
                nvariables(equations) * nnodes(mortar_l2)}
  fstar_upper_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
  fstar_lower_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
  u_threaded =           MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]

  (; fstar_upper_threaded, fstar_lower_threaded, u_threaded)
end

end # @muladd
