function [idx, val] = localMin(A)
%LOCALMIN Find local minimum in vector A
%   Detailed explanation goes here
minCond = [false, (A(2:end-1) < A(1:end-2)) & (A(2:end-1) < A(3:end)), false];
idx = find(minCond);
val = A(idx);
end

