
exp_group(G, p, X) = translate(G, p, exp_lie(G, translate_to_id(G, p, X, LeftSide())), (LeftAction(), LeftSide()))
exp_group!(G, tmp, p, X) = begin
    ξ = translate_to_id(G, p, X, LeftSide())
    return translate!(G, tmp, p, exp_lie(G, ξ), (LeftAction(), LeftSide()))
end
log_group(G, p, q) = translate_from_id(G, p, log_lie(G, translate(G, p, q, (RightAction(), LeftSide()))), LeftSide())
log_group!(G, tmp, p, q) = begin
    pinvq = translate(G, p, q, (RightAction(), LeftSide()))
    ξ = log_lie(G, pinvq)
    return translate_from_id!(G, tmp, p, ξ, LeftSide())
end
