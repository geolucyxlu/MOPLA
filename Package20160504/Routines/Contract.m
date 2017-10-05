function m = Contract(X,Y)
% Contract.m
% Double-index contraction between two 4th-order tensors.

% Input:  X, Y, two 4th-order tensors;

% Output: m, a 4th-order tensor;
%--------------------------------------------------------------------------

    m = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    x = reshape(X(i,j,:,:),3,3);
                    y = Y(:,:,k,l);
                    m(i,j,k,l) = sum(sum(x.*y));
                end
            end
        end
    end
end
