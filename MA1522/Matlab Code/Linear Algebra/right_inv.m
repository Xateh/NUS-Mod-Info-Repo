function [mat] = right_inv(mat)
%Right_INV Summary of this function goes here
%   Detailed explanation goes here
[r, c] = size(mat);
if rank(mat) == r && c >= r
    mat = mat' * (mat * mat')^-1;
else
    error('No Right Inverse: Either Rank(Mat) ~= r or c < r')
end