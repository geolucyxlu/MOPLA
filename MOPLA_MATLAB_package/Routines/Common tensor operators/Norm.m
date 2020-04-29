function m = Norm(X)
% Norm.m
% The Norm of a 4th order tensor

% Input: X: a 4th-order tensor, 3*3*3*3 matrix;

% Output: a scalar,indicating the magnitude of the 4th order tensor.
    m=0;
    for i=1:3
        for j=1:3
            m = m+norm(reshape(X(i,j,:,:),3,3),'fro');
        end
    end
end