function [mat] = ero(instruct, size, row1, row2, amount)
%ERO Summary of this function goes here
%   Detailed explanation goes here

arguments
   instruct;
   size;
   row1;
   row2 = row1;
   amount = 1;
end

mat = eye(size);
switch instruct
    case 'swap'
        mat([row1, row2],:) = mat([row2, row1],:);
    case 'add'
        mat(row1,:) = mat(row1,:) + amount*mat(row2,:);
    case 'multiply'
        mat(row1,:) = amount*mat(row1,:);
    otherwise
        error('Invalid Input');
end