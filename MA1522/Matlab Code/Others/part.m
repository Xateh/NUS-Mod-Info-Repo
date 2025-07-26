% Define the symbolic variable
syms x y;

% Define the rational function
% R = (3*x - 6)/(4 - 2*x);
% R = 3/x;
R = 1/(x^6 + 1);

% Perform partial fraction decomposition
R = simplify(R);
decomp = partfrac(R)
% num_terms = numel(regexp(char(decomp), '[+-]?'))+1;
num_terms = numel(decomp)

% Extract terms
if num_terms == 1 && decomp == R
    terms = {R};
else
    terms = children(decomp);
end
disp(terms);
% Initialize containers for polynomial and fractional terms
polynomial_terms = [];
fractional_terms = [];

function bool = is_zero(expr)
    symbolicVars = symvar(expr);
    if isempty(symbolicVars)
        if expr == 0
            bool = true;
        else
            bool = false;
        end
        % bool = expr == 0;
    else
        sol = solve(expr == 0, symbolicVars);
        % if isempty(sol)
        %     bool = false;
        % else
        %     bool = true;
        % end
        bool = ~isempty(sol);
    end
end

% [num, den] = numden(sym(1/x));
% is_zero(sym(1))


% Analyze each term
for i = 1:length(terms)
    term = terms{i};

    % Extract numerator and denominator
    [num, den] = numden(term);
    fprintf('%s     %s \n', string(num), string(den));
    if is_zero(den)
        fractional_terms = [fractional_terms; term];
    else
        polynomial_terms = [polynomial_terms; term];
    end
end

% Display polynomial terms
disp('Polynomial terms:');
disp(polynomial_terms);

% Display fractional terms
disp('Fractional terms:');
disp(fractional_terms);
