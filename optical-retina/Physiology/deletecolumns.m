function [slimMatrix] = deletecolumns (matrixName, columnsToBeDeleted)
% This function deletes the indexed columns from a matrix. Enter
% columnsToBeDeleted as a [row vector] (e.g [1 4 5] or [2:6]).

[m n] = size (matrixName);
keepVector(1:n) = 1;
keepVector(columnsToBeDeleted) = 0;
keepMatrix = diag(keepVector);
slimMatrix = matrixName * keepMatrix;
slimMatrix(:,~any(slimMatrix)) = [];