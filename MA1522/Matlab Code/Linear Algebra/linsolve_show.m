
function [output, rref_mat, x, row_ops, elem_matrices] = linsolve_show(A, B, divide_zero, show_steps, show_steps_mat)
% Input: A is the augmented matrix, B is the B of Matrix Equation Ax = B
% Output: rref_mat is the RREF of A
% row_ops is a cell array containing the row operations performed to obtain the RREF form
% elem_matrices is a cell array containing the corresponding elementary matrices for each row operation
% Initialization
arguments
    A;
    B;
    divide_zero = false;
    show_steps = false;
    show_steps_mat = show_steps;
end
[ra, ca] = size(A);
[rb, cb] = size(B);
output = [];
rref_mat = [];
x = [];
row_ops = {};
elem_matrices = {};

% Check if able to solve by matrix dimension observation
if ra ~= rb
    error('No. of Rows of A (%d) not equal to No. of Rows of B (%d)', ra, rb);
end

% Get information from rref_show(A) function
[rref_mat, row_ops, elem_matrices] = rref_show(A, divide_zero,false, false);

% Get x by transforming B
x = B;
step_mat = A;
for i = 1:length(row_ops)
    x = elem_matrices{i} * x;
    % Print the row operations and corresponding elementary matrices
    if show_steps
        fprintf('Step %d: %s\n',i,row_ops{i});
    end

    if show_steps_mat
        step_mat = elem_matrices{i} * step_mat;
        disp([step_mat x]);
    end
end

x = num2cell(x);
[rx, cx] = size(x);
x = reshape([x{:}], rx, cx);
output = [rref_mat, x];

end
% Compute RREF and row operations