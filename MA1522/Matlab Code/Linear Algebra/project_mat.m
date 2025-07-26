function [projection] = project_mat(subspace, vector)
%PROJECT_MAT Summary of this function goes here
%   Detailed explanation goes here

projection = zeros(size(vector, 1), 1);
orthbasis = orth(subspace);
for c = 1:size(orthbasis, 2)
    colvector = orthbasis(:, c);
    projection = projection + dot(vector, colvector)/norm(colvector)^2*colvector;
end

% alternatively, find lsqr(subspace, vector) then premultiply with subspace
% projection = subspace * lsqr(subspace, vector);

end