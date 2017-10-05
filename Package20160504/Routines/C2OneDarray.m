function C1d = C2OneDarray(C)
% C2OneDarray.m
% convert 4th-order C(i,j,k,l) to 1d array(1x81) row-major
% i,j,k,l: 1-3
% 
% Input:  C, 3*3*3*3 matrix;
%
% Output: C1d, 1*81 matrix.
%--------------------------------------------------------------------------
C1d = zeros(1,81);
jj = 1;

for k=1:3
   for l=1:3
      for i=1:3
         for j=1:3
             C1d(1,jj)=C(i,j,k,l);
             jj = 1+jj;
         end
      end
   end
end
end