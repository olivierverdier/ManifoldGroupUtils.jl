module ManifoldGroupUtils

import ManifoldsBase
using Manifolds
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

include("rotation_action.jl")


end
