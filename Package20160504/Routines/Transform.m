function ss=Transform(X,Q)
% Transform.m    
% 4th-order tensor transformation between coordinate systems. 
%
% Input:  X, 3*3*3*3 matrix;
%         Q, 3*3 matrix;

% Output: ss, 3*3*3*3 matrix.
%--------------------------------------------------------------------------
    ss = zeros(3,3,3,3);
    qq = zeros(3,3,3,3);
    
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    for m=1:3
                        for n=1:3
                            for s=1:3
                                for t=1:3
                                    qq(m,n,s,t)= Q(i,m)*Q(j,n)*Q(k,s)*Q(l,t);
                                end
                            end
                        end
                    end
                    ss(i,j,k,l)=sum(reshape(qq,[],1).*reshape(X,[],1));
                end
            end
        end
    end
end