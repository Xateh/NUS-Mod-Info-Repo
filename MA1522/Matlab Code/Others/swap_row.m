function A = swap_row(A, i, j, verbosity)
    arguments
        A 
        i 
        j 
        verbosity = 0
    end
    A([j, i], :) = A([i, j], :);
    A = simplify(A);
    if verbosity >= 1
        fprintf('R%d <-> R%d \n', i, j);
    end
    if verbosity >= 2
        disp(A);
    end
end