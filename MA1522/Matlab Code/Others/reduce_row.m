function A = reduce_row(A, i, scalar, j, verbosity)
    arguments
        A 
        i 
        scalar 
        j 
        verbosity = 0
    end
    [~, den] = numden(scalar);
    if is_zero(den)
        fprintf('Error: Attempted to scale Row%d by Inf \n', i);
    end

    A(i, :) = A(i, :) - scalar * A(j, :);
    A = simplify(A);
    if verbosity >= 1
        fprintf('R%d -> R%d - (%s)R%d \n', i, i, string(scalar), j);
    end
    if verbosity >= 2
        disp(A);
    end
end