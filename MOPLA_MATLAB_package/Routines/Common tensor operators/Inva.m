function m = Inva(x)
% Inva.m
% The invariant of a 2nd order tensor

% Input: x: a 2nd-order tensor, 3*3 matrix;

% Output: m: the invariant, a scalar

    m = (0.5* contract1(x,x))^0.5;
end

function m = contract1(x,y)
    m = sum(sum(x.*y));
end