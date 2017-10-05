function m = contract1(x,y)
% contract1.m
% Double-index contraction between two 2nd-order tensors.
%--------------------------------------------------------------------------
    m = sum(sum(x.*y));
end