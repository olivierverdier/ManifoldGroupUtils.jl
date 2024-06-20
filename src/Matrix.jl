
"""
    get_id_matrix_lie(G)

The identity matrix on the Lie algebra of the group `G`.
"""
function get_id_matrix_lie(G) 
    dim = manifold_dimension(G)
    T = ManifoldsBase.allocate_result_type(G, typeof(get_id_matrix_lie), ())
    return Matrix{T}(LinearAlgebra.I, dim, dim)
end


"""
    compose_matrix_op(
      G, # group
      M,p, # Manifold, point
      op, # operator Alg(G) -> T_pM
      mat, # matrix in basis BG
      BG, # basis of Alg(G)
      BM, # basis of T_pM
      )

Compose a matrix `mat` of a linear endomorphism of Alg(G)
in some basis `BG` with an operator op : ``Alg(G) → T_p M``
itself equipped with a basis `BM`.
"""
function compose_matrix_op(
    G, # group
    M,p, # Manifold, point
    op, # operator Alg(G) -> T_pM
    mat, # matrix in basis BG
    BG::AbstractBasis, # basis of Alg(G)
    BM::AbstractBasis, # basis of T_pM
    )
    idim = size(mat,2)
    odim = manifold_dimension(M)
    T = ManifoldsBase.allocate_result_type(G, typeof(compose_matrix_op), ())
    rmat = Array{T}(undef, odim, idim)
    for (b, rb) in zip(eachcol(mat), eachcol(rmat))
        vec = op(get_vector_lie(G, b, BG))
        rb[:] = get_coordinates(M, p, vec, BM)
    end
    return rmat
end

"""
    get_op_matrix(G, # group
      M, p, # manifold + point
      op, # operator from Alg(G) -> T_pM
      BG, # Lie algebra basis
      BM, # basis at T_pM
      )

Matrix of an operator op : ``Alg(G) → T_p M``
computed in the basis `BG` and `BM`.
"""
get_op_matrix(G, # group
                        M, p, # manifold + point
                        op, # operator from Alg(G) -> T_pM
                        BG::AbstractBasis, # Lie algebra basis
                        BM::AbstractBasis, # basis at T_pM
                        ) = compose_matrix_op(G, M, p, op, get_id_matrix_lie(G), BG, BM)

"""
    compose_lie_matrix_op(
        G, # group
        op, # operator Alg(G) -> Alg(G)
        mat, # matrix in the basis
        B # basis of Alg(G)
        )

Compute the matrix in the basis `B`
of the composition of `op` a linear endomorphism  of Alg(G)
and a matrix `mat`, also expressed in the basis `B`.
"""
compose_lie_matrix_op(
    G, # group
    op, # operator Alg(G) -> Alg(G)
    mat, # matrix in the basis
    B::AbstractBasis # basis of Alg(G)
    ) = compose_matrix_op(G, G, identity_element(G), op, mat, B, B)


"""
    matrix_from_lin_endomorphism(
        G, # group
        op, # Alg(G) -> Alg(G)
        B # basis of Alg(G)
        )

Compute the matrix in the basis `B`
 of an operator `op` on a Lie algebra Alg(G).
"""
matrix_from_lin_endomorphism(
    G, # group
    op, # Alg(G) -> Alg(G)
    B::AbstractBasis # basis of Alg(G)
    ) = compose_lie_matrix_op(G, op, get_id_matrix_lie(G), B)


"""
    get_proj_matrix(A::GroupAction, x, BG, BM)

From a group action ``G ⊂ Diff(M)``,
and a point ``x ∈ M``,
compute the projection matrix in the basis `BG` of Alg(G)
and `BM` of ``T_x M``
of the operator ``ξ ↦ ξ ⋅x``, where ``⋅`` denotes the
infinitesimal group action above.
"""
get_proj_matrix(A::AbstractGroupAction, x, BG, BM) = begin
    G = base_group(A)
    op(ξ) = apply_diff_group(A, Identity(G), ξ, x)
    mat = GroupTools.get_id_matrix_lie(G)
    return GroupTools.compose_matrix_op(G, group_manifold(A), x, op, mat, BG, BM)
end
