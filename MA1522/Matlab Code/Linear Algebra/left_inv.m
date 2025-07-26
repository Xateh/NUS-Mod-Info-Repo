function [mat] = left_inv(mat)
%LEFT_INV Summary of this function goes here
%   Detailed explanation goes here
[r, c] = size(mat);
if rank(mat) == c && r >= c
    mat = ((mat' * mat)^-1) * mat';
else
    error('No Left Inverse: Either Rank(Mat) ~= c or r < c')
end