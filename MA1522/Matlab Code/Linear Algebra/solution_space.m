function [solution_cell, solution, solution_space] = solution_space(A, b)
% Generate solution for the equation A*x == b

%{
Detailed explanation goes here

solution_space returns cell arrays specifying the solution space for A*x == b with conditions.

Input Argument:
A (necessary)
b (optional): defaulted to zero vector.

Output:
solution_cell: cell array of all possible solutions
solution: struct of all possible solutions
solution_space (NOT IMPLEMENTED YET): cell array of possible solutions by space basis
%}

arguments
    A;
    b = sym(zeros(size(A, 1),1));
end

% Display A, b for visual checks of Input Arguments
display(A);
display(b);

% Initialise required variables
solution_space = cell(0);
solution_mat = sym('x', [size(A, 2) 1]);
symlist = union(solution_mat, [symvar(sym(A)), symvar(sym(b))], 'stable');

% Solve
solution = solve(A*solution_mat == b, ...
                symlist, ...
                ReturnConditions=true);

solution_cell = solution_struct2cell(solution);

end

% Custom struct2cell for struct output of solve function
function [solution_cell] = solution_struct2cell(sol)
    fn = fieldnames(sol);
    solution_cell = cell(numel(fn),2);
    
    for fieldnum = 1:numel(fn)
        field = sol.(fn{fieldnum});
        
        solution_cell{fieldnum,1} = fn{fieldnum};
        if ~isempty(field)
            solution_cell(fieldnum,2:1+numel(field)) = sym2cell(field);
        end
    end
end