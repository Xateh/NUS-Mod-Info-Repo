
function [rref_mat, row_ops, elem_matrices] = rref_show(A, divide_zero, show_steps, show_steps_mat)
% Input: A is the augmented matrix
% divide_zero (default: false): whether to divide variables that could be
% zero
% show_steps (default: true): whether to display steps
% show_steps_mat (default: show_steps): whether to display matrix at
% current state
% Output: rref_mat is the RREF of A
% row_ops is a cell array containing the row operations performed to obtain the RREF form
% elem_matrices is a cell array containing the corresponding elementary matrices for each row operation
% Initialization
arguments
    A;
    divide_zero = false;
    show_steps = true;
    show_steps_mat = show_steps;
end

[r, c] = size(A);
rref_mat = A;
row_ops = cell(r * r, 1);
step_mat = cell(r * r, 1);
elem_matrices = cell(r * r, 1);
disp(A);
% Perform Gaussian elimination to obtain upper triangular matrix
for i = 1:min(r,c)
    passthrough_mat = rref_passthrough(rref_mat, i);
    if isa(rref_mat,'sym')
        while ~all(all(isAlways(rref_mat == passthrough_mat, "Unknown", false))) && ~isa(passthrough_mat, 'char')
            rref_mat = passthrough_mat;
            passthrough_mat = rref_passthrough(rref_mat, i);
        end
    elseif ~isa(passthrough_mat,'char')
        rref_mat = passthrough_mat;
    end
end

% Perform back-substitution to obtain reduced row echelon form
for i = min(r,c):-1:2
    if rref_mat(i,i) == 0
        continue
    end
    for j = i-1:-1:1
        if rref_mat(j, i) ~= 0 && row_pivot_column(rref_mat(i,:)) == i
            rref_mat = row_elimination_and_record(rref_mat, j, i);
        end
    end
end

% Shift Rows by Column
for i = 1:r
    rpc = row_pivot_column(rref_mat(i,:));
    curr = i;
    while curr ~= 1
        if rpc < row_pivot_column(rref_mat(curr - 1,:))
            rref_mat = row_swap_and_record(rref_mat, curr, curr-1);
            curr = curr - 1;
        else
            break;
        end
    end
end

% Round the entries to 6 decimal places to avoid rounding errors
% rref_mat = eval(rref_mat);
% rref_mat = double(rref_mat);
if ~isa(rref_mat, 'sym')
    rref_mat = double(rref_mat);
    rref_mat = round(rref_mat, 6);
end

% Remove Empty Cells
row_ops = row_ops(~cellfun('isempty', row_ops));
elem_matrices = elem_matrices(~cellfun('isempty', elem_matrices));

% Print the row operations and corresponding elementary matrices
%{
disp(A);
for i = 1:length(row_ops)
    if show_steps
        fprintf('Step %d: %s\n',i,row_ops{i});
    end

    if show_steps_mat
        disp(step_mat{i});
    end
end
%}

function [mat] = rref_passthrough(rref_mat, i)
    mat = rref_mat;
    % Find the first non-zero entry in the ith column
    [pivot_row_value, pivot_row] = select_pivot_row(i, i);

    % If there is no non-zero entry in the ith column, move to the next column
    if pivot_row_value == 0
        mat = 'next';
        return;
    end
    % Swap rows if necessary to bring pivot to the ith row
    if pivot_row > i
        mat = row_swap_and_record(mat, pivot_row, i);
    end
    % Scale the ith row to make the pivot equal to 1
    if mat(i, i) ~= 1
        mat = row_unit_and_record(mat, i);
    end
    % Eliminate the entries below the pivot in the ith column
    for j = i+1:r
        if mat(j, i) ~= 0
            mat = row_elimination_and_record(mat, j, i);
        end
    end
end

function [pivot_row_value, pivot_row] = select_pivot_row(start_row, col)
    pivot_row = 0;
    pivot_row_value = sym(0);
    
    for row = start_row:r
        row_value = rref_mat(row, col);
        if ~isAlways(row_value == 0, 'Unknown', false)
            if isa(pivot_row_value, 'sym') && isa(row_value, 'sym') && ~isAlways(pivot_row_value == row_value, "Unknown", false)
                if (has(pivot_row_value, row_value) && hasSymType(pivot_row_value, 'variable') && hasSymType(row_value, 'variable')) || ...
                        (~isAlways(pivot_row_value == 1, "Unknown", false) && isAlways(row_value == 1, "Unknown", false))
                    pivot_row = row;
                    pivot_row_value = row_value;
                end
            end
            if isAlways(pivot_row_value == 0, 'Unknown', false)
                pivot_row = row;
                pivot_row_value = row_value;
            elseif isAlways(pivot_row_value > row_value, 'Unknown', false) && ~isAlways(pivot_row_value == 1, "Unknown", false)
                pivot_row = row;
                pivot_row_value = row_value;
            elseif isAlways(pivot_row_value < 0, 'Unknown', false) && isAlways(row_value > 0, "Unknown", false)
                pivot_row = row;
                pivot_row_value = row_value;
            end
        end
    end
end

function [rpc] = row_pivot_column(row)
    rpc = 1;
    done = false;
    for index = 1:length(row)
        item = row(index);
        if isa(item, 'sym')
            if ~isequaln(item, sym(0))
                done = true;
            end
        else
            if item ~= 0
                done = true;
            end
        end
        if done
            break;
        end
        rpc = rpc + 1;
    end
end

%{
function [bool] = is_zero_row(row)
    bool = true;
    for index = 1:length(row)
        if row(index) ~= 0
            bool = false;
        end
    end
end
%}

function [mat] = row_swap_and_record(mat, row_1, row_2)
    step_num = nnz(~cellfun(@isempty, step_mat)) + 1;
    mat([row_2 row_1], :) = mat([row_1, row_2], :);
    if isa(mat,'sym'); mat = simplify(mat); end
    step_mat{step_num} = mat;
    next_step = sprintf('R%d <-> R%d', row_1, row_2);
    row_ops{step_num} = next_step;
    elem_matrice_to_add = sym(eye(r));
    elem_matrice_to_add([row_2 row_1], :) = elem_matrice_to_add([row_1 row_2], :);
    elem_matrices{step_num} = elem_matrice_to_add;

    if show_steps; fprintf('Step %d: %s\n', step_num, next_step); end
    if show_steps_mat; disp(mat); end
end

function [mat] = row_unit_and_record(mat, iv)
    pivot = mat(iv, iv);
    if isa(mat,"sym") && ~divide_zero
        if hasSymType(pivot, "variable")
            coeff_list = coeffs(pivot);
            if coeff_list(end) >= 0
                return;
            end
            pivot = -1;
        end
    end
    step_num = nnz(~cellfun(@isempty, step_mat)) + 1;
    mat(iv, :) = mat(iv, :) / pivot;
    if isa(mat,'sym'); mat = simplify(mat); end
    step_mat{step_num} = mat;
    next_step = sprintf('R%d -> (1/%s) R%d', iv, string(pivot), iv);
    row_ops{step_num} = next_step;
    elem_matrice_to_add = sym(eye(r));
    elem_matrice_to_add(iv, iv) = 1/pivot;
    elem_matrices{step_num} = elem_matrice_to_add;

    if show_steps; fprintf('Step %d: %s\n', step_num, next_step); end
    if show_steps_mat; disp(mat); end
end

function [mat] = row_elimination_and_record(mat, jv, iv)
    step_num = nnz(~cellfun(@isempty, step_mat)) + 1;
    factor = 0;
    if mat(iv, iv) == 1
        factor = mat(jv, iv);
    elseif isa(mat,'sym')
        if has(mat(jv, iv), mat(iv, iv))
            subeq = subs(mat(jv,iv),mat(iv,iv),sym('substitute_value'));
            coeff_list = coeffs(subeq + 1, sym('substitute_value'));
            for coeff = 1:length(coeff_list)
                if coeff > 1
                    factor = factor + coeff_list(coeff)*mat(iv, iv)^(coeff - 2);
                end
            end
        end
    end

    if isAlways(factor == 0, "Unknown", false)
        return
    end

    mat(jv, :) = mat(jv, :) - factor*mat(iv, :);
    if isa(mat,'sym'); mat = simplify(mat); end
    step_mat{step_num} = mat;
    next_step = sprintf('R%d -> R%d - (%s) R%d', jv, jv, string(factor), iv);
    row_ops{step_num} = next_step;
    elem_matrice_to_add = sym(eye(r));
    elem_matrice_to_add(jv, iv) = -factor;
    elem_matrices{step_num} = elem_matrice_to_add;

    if show_steps; fprintf('Step %d: %s\n', step_num, next_step); end
    if show_steps_mat; disp(mat); end
end

end
% Compute RREF and row operations