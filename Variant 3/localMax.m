function [idx, val] = localMax(A)
%LOCALMAX Find local maximum in vector A
%   Detailed explanation goes here
maxCond = [false, (A(2:end-1) > A(1:end-2)) & (A(2:end-1) > A(3:end)), false];
idx = find(maxCond);
val = A(idx);
end

