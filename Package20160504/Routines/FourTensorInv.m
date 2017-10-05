function m = FourTensorInv(X)
% FourTensorInv.m
% Inverse of a 4th-order symmetric Tensor 
%--------------------------------------------------------------------------
    
  sqrt2 = 2^0.5;
  sqrt3 = 3^0.5;
  
b(:,:,1) = 1/(sqrt2*sqrt3) * [-1, 0, 0; 0, -1, 0; 0, 0, 2];
    b(:,:,2) = 1/sqrt2 * [-1, 0, 0; 0, 1, 0; 0, 0, 0];
    b(:,:,3) = 1/sqrt2 * [0, 0, 0; 0, 0, 1; 0, 1, 0];
    b(:,:,4) = 1/sqrt2 * [0, 0, 1; 0, 0, 0; 1, 0, 0];
    b(:,:,5) = 1/sqrt2 * [0, 1, 0; 1, 0, 0; 0, 0, 0];
    b(:,:,6) = 1/sqrt3 * [1, 0, 0; 0, 1, 0; 0, 0, 1]; 
    
    h = 5;
    M = zeros(h,h);
    for lambda = 1:h
        for xi = 1:h
            M (lambda,xi)= contract1(Multiply(X,b(:,:,xi)),b(:,:,lambda));  
        end
    end
    
    m = ConvertBK(eye(5)/M,h,b);
end

function M = ConvertBK(x,h,b)
    M = zeros(3,3,3,3);
    bb = zeros(h,h);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                     for m=1:h
                         for n=1:h
                             bb(m,n)= b(i,j,m)*b(k,l,n); 
                         end
                     end
                   M(i,j,k,l)=sum(sum(x.*bb));
                end
            end
        end
    end
end
