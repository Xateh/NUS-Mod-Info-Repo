function bool = is_scalable(expr)
    [num, den] = numden(expr);
    bool = ~is_zero(num) && ~is_zero(den);
end
