
translate_diff_id(G, χ, ξ, conv) = translate_diff(G, χ, identity_element(G), ξ, conv)

"""
The left translation ``T_L(g,ξ) = gξ``.
"""
move(G, χ, ξ, ::LeftSide) = translate_diff_id(G, χ, ξ, (LeftAction(), LeftSide()))
"""
The right translation ``T_R(g,ξ) = ξg``.
"""
move(G, χ, ξ, ::RightSide) = translate_diff_id(G, χ, ξ, (RightAction(), RightSide()))

#--------------------------------

_get_side(::LeftAction) = RightSide()
_get_side(::RightAction) = LeftSide()


_transporter(G, χ, ξ, dir) = move(G, χ, ξ, _get_side(dir))

"""
This should hold for *any* group action ``A`` on any manifold.
If you define ``π_x(g) := A(g, x)`` for ``g ∈ G`` and ``x ∈ M``,
and define, for ``ξ in Alg(G)``,
 ``T_R(g, ξ) := ξg`` (the right translation),
and ``T_L(g, ξ) := gξ`` (the left translation), then we have the identity:
```math
⟨Dπ_{x}(g), T(g, ξ)⟩ = ⟨Dπ_{A(g,x)}(1), ξ⟩
```
where, for a *left* action, ``T`` is the *right* translation,
and for a *right* action, ``T`` is the *left* translation.
"""
check_apply_diff_group(A::AbstractGroupAction{TAD}, χ, ξ, p) where {TAD} = begin
    G = base_group(A)
    p_ = apply(A, χ, p)
    v1 = apply_diff_group(A, χ, _transporter(G, χ, ξ, TAD()), p)
    v2 = apply_diff_group(A, identity_element(G), ξ, p_)
    return isapprox(TangentSpace(G, p_), v1, v2)
end

#--------------------------------

"""
Test the differential of the inverse on a Lie group.
Denote this inverse by ``I(g) := g^{-1}``.
If the left and right transports are ``T_L(g,ξ) := gξ``
and ``T_R(g,ξ) := ξg`` respectively, then
```math
⟨DI(g), T_L(g,ξ)⟩ = -T_R(g^{-1}, ξ)
```
and
``` math
⟨DI(g), T_R(g,ξ)⟩ = -T_L(g^{-1}, ξ)
```
"""
check_inv_diff(G, χ, ξ, conv) = begin
    χ_ = inv(G, χ)
    computed = inv_diff(G, χ, move(G, χ, ξ, conv))
    expected = -move(G, χ_, ξ, switch_side(conv))
    return isapprox(TangentSpace(G, χ_), computed, expected)
end
