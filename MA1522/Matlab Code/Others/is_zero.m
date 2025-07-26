function bool = is_zero(expr)
    symbolicVars = symvar(expr);
    if isempty(symbolicVars)
        if expr == 0
            bool = true;
        else
            bool = false;
        end
    else
        sol = solve(expr == 0, symbolicVars);
        bool = ~isempty(sol);
    end
end