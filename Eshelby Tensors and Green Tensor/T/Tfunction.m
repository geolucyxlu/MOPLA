function T = Tfunction(x, Cm, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3, Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww, type)
% Tfunction.m
% 4th-order Green tensor T
% Qu et al., 2016 
% Input:  x, three semi-axes of the inclusion, 3*1 matrix;
%         Cm, stiffness of the matrix, 4th order matrix;
%         Alp#, Bet#, ww#, nodes and weights;
%         type,           1 for compressible material
%                         2 for incompressible material
%
% Output: T, 4th-order tensor, 3*3*3*3 matrix;
%% updates:
%  30th July, 2019 by Lucy: Extended to compressible material 
%                           Combine the C2OneDarray function  
%--------------------------------------------------------------------------
%  rewrite the 4th order matrix stiffness tensor into a 1D array format
   C = C2OneDarray(Cm);
%  get PHI and R according to three semi-axes of the inclusion 
   [PHI, R]      = Incl(x);
%  compute r according to Empirical function.  Qu et al., 2016
   [r0,r1,r2,r3] = EmprFun(PHI);
%  allocate T before the loop   
   T             = zeros(3,3,3,3);
   
%  according to the estimated CPUtime, choose the optimal quadrature for T  
   for i = 1:3
       for j = 1:3
           for k = i:3
               for l = j:3
                    if (R<r0)


                          if PHI < 45				 
                             T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp2, Bet2, ww2, C, type); 
                          else
                             T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp1, Bet1, ww1, C, type); 
                          end

                    elseif (R >= r0) && (R < r1)

                          if PHI < 15
                             T(i,j,k,l) = GGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, Alp4,...
                                             Bet4, ww4, C, type);
                          else
                             T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp3, Bet3, ww3, C, type); 
                          end

                    elseif (R >= r1) && (R < r2)

                          if PHI < 65
                             T(i,j,k,l) = GGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, Alp5,...
                                             Bet5, ww5, C, type);
                          else
                             T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp3, Bet3, ww3, C,type);
                          end

                    elseif (R >= r2) && (R < r3)

                          if PHI < 70
                             T(i,j,k,l) = GGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, Alp6,...
                                             Bet6, ww6, C, type);
                          else
                             T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp3, Bet3, ww3, C, type);
                          end

                    else   
                           T(i,j,k,l) = AGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, p,...
                                             ww, 10^(-4), C, 1, type);
                    end

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