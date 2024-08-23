module ManifoldGroupUtils

import ManifoldsBase
using Manifolds
import LinearAlgebra # just for `I`, the identity object

include("Matrix.jl")

export algebra, rand_lie, translate_to_id, translate_from_id

"""
    algebra(G)

The tangent space at identity of the group G.
"""
algebra(G) = TangentSpace(G, identity_element(G))

#$ TODO: should be adjoint_action(G, p, X, RightAction())
inverse_adjoint_action(G::AbstractDecoratorManifold, p, X) = adjoint_action(G, inv(G, p), X)

translate_diff_id(G, χ, ξ, conv) = translate_diff(G, χ, identity_element(G), ξ, conv)

"""
    translate_to_id(G, χ, v, ::GroupActionSide)

Compute ``η = v χ⁻¹``, or ``χ⁻¹ v`` depending on whether the group action side is `Right` or `Left` respectively.
"""
translate_to_id(G, χ, v, ::RightSide) = inverse_translate_diff(G, χ, χ, v, (RightAction(), RightSide()))
translate_to_id(G, χ, v, ::LeftSide) = inverse_translate_diff(G, χ, χ, v, (LeftAction(), LeftSide()))



"""
    translate_from_id(G, χ, ξ, ::GroupActionSide)

The left translation ``T_L(g,ξ) = gξ`` for the `Left` side,
 or right translation ``T_R(g,ξ) = ξg`` for the `Right` side.
"""
translate_from_id(G, χ, ξ, ::LeftSide) = translate_diff_id(G, χ, ξ, (LeftAction(), LeftSide()))
translate_from_id(G, χ, ξ, ::RightSide) = translate_diff_id(G, χ, ξ, (RightAction(), RightSide()))


"""
    rand_lie(rng::AbstractRNG, G)

Random element in the Lie algebra of the group `G`.
"""
rand_lie(rng, G) = rand(rng, algebra(G))

include("rotation_action.jl")


end
