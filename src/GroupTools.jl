module GroupTools

import ManifoldsBase
import Manifolds:
    TangentSpace,
    manifold_dimension,
    AbstractBasis, get_vector_lie, get_coordinates,
    identity_element
import Manifolds:
    AbstractGroupAction, base_group, group_manifold,
    LeftAction, RightAction, LeftSide, RightSide, apply, apply_diff_group,
    translate_diff, inv_diff, switch_side
import LinearAlgebra # just for `I`, the identity object

include("Matrix.jl")
include("TestUtils.jl")

algebra(G) = TangentSpace(G, identity_element(G))

end
