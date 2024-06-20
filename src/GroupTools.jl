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
include("TestUtils.jl")

algebra(G) = TangentSpace(G, identity_element(G))
inverse_adjoint_action(G::AbstractDecoratorManifold, p, X) = adjoint_action(G, inv(G, p), X)

end
