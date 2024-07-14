using ManifoldGroupUtils
using Test
using Manifolds
import ManifoldsBase

import Random
rng = Random.default_rng()


compose_matrix_op_(G,M,p,op,mat,BG,BM) = begin
    vecs = [op(get_vector_lie(G, b, BG)) for b in eachcol(mat)]
    cols = [get_coordinates(M, p, vec, BM) for vec in vecs]
    return hcat(cols...)
end

check_compose_matrix_op(A, mat, p) = begin
    G = base_group(A)
    M = group_manifold(A)
    op(ξ) = apply_diff_group(A, Identity(G), ξ, p)
    args = G, M, p, op, mat, DefaultOrthogonalBasis(), DefaultOrthogonalBasis()
    computed = ManifoldGroupUtils.compose_matrix_op(args...)
    expected = compose_matrix_op_(args...)
    return computed ≈ expected
end

check_compose_matrix_op_rng(rng, A) = begin
    d = manifold_dimension(base_group(A))
    mat = randn(rng, d, d)
    p = rand(rng, group_manifold(A))
    return check_compose_matrix_op(A, mat, p)
end

@testset "compose_matrix_op" begin
    A = RotationAction(Euclidean(4), SpecialOrthogonal(4))
    @test check_compose_matrix_op_rng(rng, A)
end


function get_proj_matrix_(A::AbstractGroupAction, x, BG, BM)
    G = base_group(A)
    M = group_manifold(A)
    idim = manifold_dimension(G)
    odim = manifold_dimension(M)
    T = ManifoldsBase.allocate_result_type(G, typeof(get_proj_matrix_), ())
    rmat = Array{T}(undef, odim, idim)

    mat = ManifoldGroupUtils.get_id_matrix_lie(G)
    for (v, rv) in zip(eachcol(mat), eachcol(rmat))
        bvec = get_vector_lie(G, v, BG)
        mvec = apply_diff_group(A, Identity(G), bvec, x)
        rv[:] = get_coordinates(M, x, mvec, BM)
    end
    return rmat
end


check_proj_matrix(A, x) = begin
    args = A, x, DefaultOrthogonalBasis(), DefaultOrthogonalBasis()
    computed = ManifoldGroupUtils.get_proj_matrix(args...)
    expected = get_proj_matrix_(args...)
    return expected ≈ computed
end

@testset "proj_matrix" begin
    A = RotationAction(Euclidean(4), SpecialOrthogonal(4))
    x = rand(rng, group_manifold(A))
    @test check_proj_matrix(A, x)
end

get_op_matrix_(G, M, p, op, BG, BM) = begin
    basis = get_basis_lie(G, BG)
    vecs = [op(b) for b in basis]
    cols = [get_coordinates(M, p, vec, BM) for vec in vecs]
    return hcat(cols...)
end


get_basis_lie(G, B::AbstractBasis) = begin
    imat = ManifoldGroupUtils.get_id_matrix_lie(G)
    basis = [get_vector_lie(G, coord, B)
             for coord in eachcol(imat)]
    return basis
end

check_get_op_matrix(A, p) = begin
    G = base_group(A)
    M = group_manifold(A)
    op(ξ) = apply_diff_group(A, Identity(G), ξ, p)
    args = G, M, p, op, DefaultOrthogonalBasis(), DefaultOrthogonalBasis()
    computed = ManifoldGroupUtils.get_op_matrix(args...)
    expected = get_op_matrix_(args...)
    return computed ≈ expected
end

@testset "get_op_matrix" begin
    A = RotationAction(Euclidean(4), SpecialOrthogonal(4))
    p = rand(rng, group_manifold(A))
    @test check_get_op_matrix(A, p)
end

compose_lie_matrix_op_(G, op, mat, B) = begin
    vecs = [op(get_vector_lie(G, c, B)) for c in eachcol(mat)]
    cols = [get_coordinates_lie(G, vec, B) for vec in vecs]
    return hcat(cols...)
end

@testset "compose_lie_matrix_op" begin
    G = SpecialOrthogonal(3)
    d = manifold_dimension(G)
    mat = rand(rng, d, d)
    mat_op = rand(rng, d, d)
    B = DefaultOrthogonalBasis()
    op(ξ) = get_vector_lie(G, mat_op*get_coordinates_lie(G, ξ, B), B)
    computed_ = compose_lie_matrix_op_(G, op, mat, B)
    computed = ManifoldGroupUtils.compose_lie_matrix_op(G, op, mat, B)
    expected = mat_op * mat
    @test computed ≈ expected
end


@testset "eltype rand_lie" begin
    G = Unitary(4)
    @test eltype(ManifoldGroupUtils.rand_lie(rng, G)) <: Complex
end
