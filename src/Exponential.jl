
exp_group(G, p, X) = apply(GroupOperationAction(G, (LeftAction(), LeftSide())), p, exp_lie(G, translate_to_id(G, p, X, LeftSide())))
exp_group!(G, tmp, p, X) = begin
    ξ = translate_to_id(G, p, X, LeftSide())
    return apply!(GroupOperationAction(G, (LeftAction(), LeftSide())), tmp, p, exp_lie(G, ξ))
end
log_group(G, p, q) = translate_from_id(G, p, log_lie(G, apply(GroupOperationAction(G, (RightAction(), LeftSide())), p, q)), LeftSide())
log_group!(G, tmp, p, q) = begin
    pinvq = apply(GroupOperationAction(G, (RightAction(), LeftSide())), p, q)
    ξ = log_lie(G, pinvq)
    return translate_from_id!(G, tmp, p, ξ, LeftSide())
end
