module ManifoldGroupUtils

import ManifoldsBase
using Manifolds
import LinearAlgebra # just for `I`, the identity object

include("Matrix.jl")

export algebra, rand_lie,
    translate_to_id, translate_to_id!,
    translate_from_id, translate_from_id!,
    exp_group, exp_group!,
    log_group, log_group!

"""
    algebra(G)

The tangent space at identity of the group G.
"""
algebra(G) = TangentSpace(G, identity_element(G))

@deprecate inverse_adjoint_action(G::AbstractDecoratorManifold, χ, X) adjoint_action(G, χ, ξ, RightAction())


"""
    translate_to_id(G, χ, v, ::GroupActionSide)

Compute ``η = v χ⁻¹``, or ``χ⁻¹ v`` depending on whether the group action side is `Right` or `Left` respectively.
"""
translate_to_id(G, χ, v, ::RightSide) = apply_diff(GroupOperationAction(G, (LeftAction(), RightSide())), χ, χ, v)
translate_to_id!(G, tmp, χ, v, ::RightSide) = apply_diff!(GroupOperationAction(G, (LeftAction(), RightSide())), tmp, χ, χ, v)
translate_to_id(G, χ, v, ::LeftSide) = apply_diff(GroupOperationAction(G, (RightAction(), LeftSide())), χ, χ, v)
translate_to_id!(G, tmp, χ, v, ::LeftSide) = apply_diff!(GroupOperationAction(G, (RightAction(), LeftSide())), tmp, χ, χ, v)



"""
    translate_from_id(G, χ, ξ, ::GroupActionSide)

The left translation ``T_L(g,ξ) = gξ`` for the `Left` side,
 or right translation ``T_R(g,ξ) = ξg`` for the `Right` side.
"""
translate_from_id(G, χ, ξ, ::LeftSide) = apply_diff(GroupOperationAction(G, (LeftAction(), LeftSide())), χ, Identity(G), ξ)
translate_from_id!(G, tmp, χ, ξ, ::LeftSide) = apply_diff!(GroupOperationAction(G, (LeftAction(), LeftSide())), tmp, χ, Identity(G), ξ)
translate_from_id(G, χ, ξ, ::RightSide) = apply_diff(GroupOperationAction(G, (RightAction(), RightSide())), χ, Identity(G), ξ)
translate_from_id!(G, tmp, χ, ξ, ::RightSide) = apply_diff!(GroupOperationAction(G, (RightAction(), RightSide())), tmp, χ, Identity(G), ξ)


"""
    rand_lie(rng::AbstractRNG, G)

Random element in the Lie algebra of the group `G`.
"""
rand_lie(rng, G) = rand(rng, algebra(G))

include("Exponential.jl")

include("rotation_action.jl")

end
