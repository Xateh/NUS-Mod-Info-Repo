function result = ref(A, verbosity)
    % Perform Gauss-Jordan elimination to reduce matrix A to REF.
    % By default, returns a matrix which is REF of A
    % If verbosity >= 1, returns {E, A} where A is the REF of A and E is
    % the product of elementary matrices such that A = E * REF_A
    arguments
        A 
        verbosity = 0;
    end

    [m, n] = size(A);
    
    if verbosity >= 1
        % Create an E matrix (to be inverted) such that A = E^-1 * REF_A
        E = sym(eye(m));
    end

    for i = 1:min(m, n) % Loop over each column
        % Find the first non-zero row in the current column
        % pivot_row = find(A(i:end, i) ~= 0, 1) + i - 1; % Find the first non-zero element in column i from row i to end
        pivot_row = get_pivot_row(A, i);

        if pivot_row == -1
        % if isempty(pivot_row)
            % If no non-zero pivot is found, continue to the next column
            continue;
        end
        
        % Swap the current row with the pivot row if necessary
        if pivot_row ~= i
            A = swap_row(A, i, pivot_row, verbosity);
            if verbosity >= 1
                E = swap_row(E, i, pivot_row, 0);
            end
        end

        % Eliminate the current column in rest of the rows below
        for j = i+1:m
            while A(j, i) ~= 0
                % simplify is necessary for comparison with decomp later on
                scalar = simplify(A(j, i) / A(i, i));
                decomp = partfrac(scalar);
                if isscalar(decomp) && decomp == scalar
                    % This runs when decomp failed to separate scalar into
                    % multiple partial fractions
                    scalars = {scalar};
                else
                    scalars = children(decomp);
                end
                for k = 1:length(scalars)
                    scalar = scalars{k};
                    [~, den] = numden(scalar);
                    % if is_scalable(scalar)
                    if ~is_zero(den)
                        A = reduce_row(A, j, scalar, i, verbosity);
                        if verbosity >= 1
                            E = reduce_row(E, j, scalar, i, 0);
                        end
                    end
                end
                scalar = simplify(A(j, i) / A(i, i));
                if A(j, i) ~= 0 && ~is_scalable(scalar)
                    A = scale_row(A, j, A(i, i), verbosity);
                    if verbosity >= 1
                        E = scale_row(A, j, A(i, i), 0);
                    end
                end
            end
        end
    end
    if verbosity >= 1
        result = {E^-1, A};
    else
        result = A;
    end
end

function pivot_row = get_pivot_row(A, col_index)
    [m, ~] = size(A);
    % Attempt to pick a pivot column that is a non-zero constant that do
    % not contain any symbols so that it is easier to reduce other rows
    for i = col_index:m
        term = A(i, col_index);
        symbolicVars = symvar(term);
        if term ~= 0 && isempty(symbolicVars)
            pivot_row = i;
            return
        end
    end
    % Attempt to pick the first non-zero row if all rows contain symbols
    for i = col_index:m
        term = A(i, col_index);
        if term ~= 0
            pivot_row = i;
            return
        end
    end
    % If the column is completely 0, return -1 to indicate non-pivot column
    pivot_row = -1;
end

% function [E, A] = swap_row(A, i, j, verbosity)
%     arguments
%         A 
%         i 
%         j 
%         verbosity = 0
%     end
%     [m, n] = size(A);
%     E = sym(eye(min(m, n)));
%     A([j, i], :) = A([i, j], :);
%     E([j, i], :) = E([j, i], :);
% 
%     if verbosity >= 1
%         fprintf('R%d <-> R%d \n', i, j);
%     end
%     if verbosity >= 2
%         disp(A);
%     end
% end

% function [E, A] = scale_row(A, i, scalar, verbosity)
%     arguments
%         A 
%         i 
%         scalar 
%         verbosity = 0
%     end
%     [m, n] = size(A);
%     E = sym(eye(m));
%     % if scalar == 0
%     %     fprintf('Error: Attempted to scale Row%d by 0 \n', i);
%     % end
%     A(i, :) = scalar * A(i, :);
%     E(i, :) = scalar * E(i, :);
% 
%     if verbosity >= 1
%         fprintf('R%d -> (%s)R%d \n', i, string(scalar), i);
%     end
%     if verbosity >= 2
%         disp(A);
%     end
% end

% function [E, A] = reduce_row(A, i, scalar, j, verbosity)
%     arguments
%         A 
%         i 
%         scalar 
%         j 
%         verbosity = 0
%     end
%     [m, n] = size(A);
%     E = sym(eye(m));
%     A(i, :) = A(i, :) - scalar * A(j, :);
%     E(i, :) = E(i, :) - scalar * E(j, :);
% 
%     if verbosity >= 1
%         fprintf('R%d -> R%d - (%s)R%d \n', i, i, string(scalar), j);
%     end
%     if verbosity >= 2
%         disp(A);
%     end
% end

% function bool = is_zero(expr)
%     symbolicVars = symvar(expr);
%     if isempty(symbolicVars)
%         if expr == 0
%             bool = true;
%         else
%             bool = false;
%         end
%     else
%         sol = solve(expr == 0, symbolicVars);
%         bool = ~isempty(sol);
%     end
% end

% syms a b c;
% AM = [a 2 a a+b a-b; 
%       a 2 a a a-b; 
%       3 3 -b 3 -b; 
%       a+1 3 a+1 a+1 a-b+1];
% get_pivot_row(AM, 1)
% result = test(AM, 2);
% A = result{1};
% A = simplify(A);
% disp(result{1});
% disp(result{2});