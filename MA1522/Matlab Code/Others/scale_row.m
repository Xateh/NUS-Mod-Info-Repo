function A = scale_row(A, i, scalar, verbosity)
    arguments
        A 
        i 
        scalar 
        verbosity = 0
    end
    [~, den] = numden(scalar);
    if scalar == 0
        fprintf('Error: Attempted to scale Row%d by 0 \n', i);
    end
    if is_zero(den)
        fprintf('Error: Attempted to scale Row%d by Inf \n', i);
    end

    A(i, :) = scalar * A(i, :);
    A = simplify(A);
    if verbosity >= 1
        fprintf('R%d -> (%s)R%d \n', i, string(scalar), i);
    end
    if verbosity >= 2
        disp(A);
    end
end