function [solution_space] = solution_space_dud(A, b)
%SOLUTION_SPACE Summary of this function goes here
%   Detailed explanation goes here

arguments
    A;
    b = 0;
end

display(A);
display(b);

solution_space = cell(0);
sym_list = union(symvar(sym(A)), symvar(sym(b)));

solution_mat = sym('x', [size(A, 2) 1]);
solution_sym = union(solution_mat, sym_list);

solution = solve(A*solution_mat == b, solution_sym, ReturnConditions=true);
sol_cell = solution_cell(solution);
display(sol_cell);

solution_space = sol_cell;

condarr = str2sym(split(string(solution.conditions), ' & '));
for condnum = 1:size(condarr,1)
    fprintf("When %s", string(~condarr(condnum)));
    % sol = solve([A*solution==b; ~condarr(condnum)], solution, ReturnConditions=true);
    new_A = subs(A, lhs(condarr(condnum)), rhs(condarr(condnum)));
    new_b = subs(b, lhs(condarr(condnum)), rhs(condarr(condnum)));
    solution = solve(new_A*solution_mat == new_b, solution_sym, ReturnConditions=true);
    
    sol_cell = solution_cell(solution);
    display(sol_cell);
    
    solution_space(:, end+1) = sol_cell(:,2:end);
end

function [display_cell] = display_solution(sol)
    fn = fieldnames(sol);
    display_cell = cell(numel(fn),2);
    
    for fieldnum = 1:numel(fn)
        field = sol.(fn{fieldnum});
        
        display_cell{fieldnum,1} = fn{fieldnum};
        if ~isempty(field)
            display_cell(fieldnum,2:1+numel(field)) = sym2cell(field);
        end
    end

    display(display_cell);
end

function [solution_cell] = solution_cell(sol)
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

function [new_condsymarr] = do_solve(A, b, condsymarr)
    new_condsymarr = cell(0);
    for condarr_row = 1:size(condarr,1)
        new_A = A;
        new_b = b;
        for condarr_col = 1:size(condarr,2)
            next_cond = condarr(condarr_row, condarr_col);
            fprintf("When %s", string(next_cond));
            % sol = solve([A*solution==b; ~condarr(condnum)], solution, ReturnConditions=true);
            if (is_equiv_cond(next_cond))
                new_A = subs(A, lhs(next_cond), rhs(next_cond));
                new_b = subs(b, lhs(condarr(condnum)), rhs(condarr(condnum)));
            else 
                new_A = [new_A; next_cond];
            end
            
        solution = solve(new_A*solution_mat == new_b, solution_mat, ReturnConditions=true);
        
        display_solution(solution);
        end
    end
end
    
function [bool] = is_no_solution(solution)
    fn = fieldnames(sol);
end

function [bool] = is_equiv_cond(cond)
    bool = contains(string(cond), '==');
end

function [conditions, cond_mat, sym_list] = cond_for_zero(A, b)
    sym_list = union(symvar(A), symvar(b));
    
    cond_mat = sym(zeros(length(sym_list) + 1, 0));
    for r = 1:size(A, 1)
        for c = 1:size(A, 2)
            if (~has(A(r,c), symvar(A)))
                continue
            end
            
            sol_zero = solve(A(r,c) == 0, sym_list, ReturnConditions=true);
            
            display(sol_zero);
            solve_mat = solve_output_matrix(sol_zero);
            
            display(solve_mat);
            
            cond_mat_T = cond_mat.';
            
            new_sol_zero_mat = solve_mat(~ismember(solve_mat, cond_mat_T, 'rows'),:);
            
            display(new_sol_zero_mat);
            if ~(isempty(new_sol_zero_mat))
                cond_mat(:, end+1) = new_sol_zero_mat.';
            end
        end
    end
    
    conditions = cell(0);
    for rnum = 1:size(cond_mat, 2)
        new_condition = sym([]);
        for symnum = 1:length(sym_list)
            symb = sym_list(symnum);
            new_condition(end+1,1) = symb == cond_mat(rnum, symnum);
        end
        
        conditions(end+1) = {new_condition};
    end
end

function [cond_mat] = solve_output_matrix(sol_struct)
    cond_mat = sym([]);
    
    f = fields(sol_struct);
    f(end-1,:) = [];
    f = string(f);
    
    for f_num = 1:length(f)
        f_name = f(f_num);
        
        cond_mat(end+1) = sol_struct.(f_name);
    end
end

end