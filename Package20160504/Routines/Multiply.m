function m = Multiply(X,y)
% Multiply.m
% Double-index contraction of a 4th-order tensor and a 2nd-order tensor.

% Input: X, a 4th-order tensor, 3*3*3*3 matrix;
%        y, a 2nd-order tensor, 3*3 matrix;

% Output: m, a 2nd-order tensor, 3*3 matrix. 
%--------------------------------------------------------------------------
    m = zeros(3,3);
    for i=1:3
        for j =1:3
             x = reshape(X(i,j,:,:),3,3);
             m(i,j) = sum(sum(x.*y));
        end
    end
end