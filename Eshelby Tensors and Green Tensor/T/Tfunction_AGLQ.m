function T = Tfunction_AGLQ(x, Cm, p, ww, type)
% TGreen.m
% 4th-order Green tensor T
% 
% Input:  x, three semi-axes of the inclusion, 3*1 matrix;
%         Cm, stiffness of the matrix, 4th order matrix;
%         p, nodes, n*1 matrix;
%         w, weights, n*1 matrix. ww = w * w';
%         type,           1 for compressible material
%                         2 for incompressible material       
% Output: T, 4th-order tensor, 3*3*3*3 matrix;
%% updates:
%  30th July, 2019 by Lucy: Extended to compressible material 
%                           Combine the C2OneDarray function  
%         type,           1 for compressible material
%                        2 for incompressible material
%--------------------------------------------------------------------------
%  rewrite the 4th order matrix stiffness tensor into a 1D array format
   C = C2OneDarray(Cm);
%  allocate T before the loop   
   T = zeros(3,3,3,3);
   
%  according to the estimated CPUtime, choose the optimal quadrature for T  
   for i = 1:3
       for j = 1:3
           for k = i:3
               for l = j:3
                    
                   T(i,j,k,l) = AGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, p,...
                                     ww, 10^(-4), C, 1, type);
                   T(k,j,i,l) = T(i,j,k,l);            
                   T(i,l,k,j) = T(i,j,k,l);
                   T(k,l,i,j) = T(i,j,k,l);                
               end
           end
       end
  end   
end

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
