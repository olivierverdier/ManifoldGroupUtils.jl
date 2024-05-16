using GroupTools
using Test
using Manifolds

@testset "GroupTools.jl" begin
    G = Orthogonal(3)
    op(ξ) = ξ
    # mat = [1. 0 0 ; 0 1 0; 0 0 1]
    mat = [1. 0  ; 0 1 ; 0 0 ]
    # TODO: test
    GroupTools.compose_matrix_op(G, G, identity_element(G), op, mat, DefaultOrthonormalBasis(), DefaultOrthonormalBasis())
end
