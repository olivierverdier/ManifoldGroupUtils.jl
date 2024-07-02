module GroupTools

import ManifoldsBase
import Manifolds:
    AbstractDecoratorManifold, TangentSpace,
    manifold_dimension,
    AbstractBasis, get_vector_lie, get_coordinates,
    identity_element
import Manifolds:
    Identity
import Manifolds:
    AbstractGroupAction, base_group, group_manifold,
    LeftAction, RightAction, LeftSide, RightSide, apply, apply_diff_group,
    translate_diff, inv_diff, switch_side,
    adjoint_action
import LinearAlgebra # just for `I`, the identity object

include("Matrix.jl")

"""
    algebra(G)

The tangent space at identity of the group G.
"""
algebra(G) = TangentSpace(G, identity_element(G))

inverse_adjoint_action(G::AbstractDecoratorManifold, p, X) = adjoint_action(G, inv(G, p), X)


"""
    rand_lie(rng::AbstractRNG, G)

Random element in the Lie algebra of the group `G`.
"""
rand_lie(rng, G) = rand(rng, algebra(G))



end
