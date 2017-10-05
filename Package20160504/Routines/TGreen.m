function T = TGreen(x, C, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3, Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww)
% TGreen.m
% 4th-order Green tensor T
% 
% Input:  x, three semi-axes of the inclusion, 3*1 matrix;
%         C, stiffness of the matrix, 1*81 matrix;
%         Alp#, Bet#, ww#, nodes and weights;
%       
% Output: T, 4th-order tensor, 3*3*3*3 matrix;
%--------------------------------------------------------------------------

%  get PHI and R according to three semi-axes of the inclusion 
   [PHI, R]      = Incl(x);
%  compute r according to Empirical function   
   [r0,r1,r2,r3] = EmprFun(PHI);
%  allocate T before the loop   
   T             = zeros(3,3,3,3);
   
%  according to the estimated CPUtime, choose the optimal quadrature for T  
   for i = 1:3
       for j = 1:3
           for k = i:3
               for l = j:3
			        if(R<r0)
                      if PHI < 45				 
				         T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp2, Bet2, ww2, C); 
                      else
					     T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp1, Bet1, ww1, C); 
                      end
                    elseif (R >= r0) && (R < r1)
                      if PHI < 15
                         T(i,j,k,l) = GGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, Alp4,...
                                         Bet4, ww4, C);
                      else
                         T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp3, Bet3, ww3, C); 
                      end
                    elseif (R >= r1) && (R < r2)
                      if PHI < 65
                         T(i,j,k,l) = GGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, Alp5,...
                                         Bet5, ww5, C);
                      else
                         T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp3, Bet3, ww3, C);
                      end
                    elseif (R >= r2) && (R < r3)
                      if PHI < 70
                         T(i,j,k,l) = GGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, Alp6,...
                                         Bet6, ww6, C);
                      else
                         T(i,j,k,l) = GLeQ([i,j,k,l], x, Alp3, Bet3, ww3, C);
                      end
                    else   
                       T(i,j,k,l) = AGLQ([i,j,k,l], x, 0, 2*pi, 0, pi, p,...
                                         ww, 10^(-4), C, 1);
                    end
                       T(k,j,i,l) = T(i,j,k,l);            
                       T(i,l,k,j) = T(i,j,k,l);
                       T(k,l,i,j) = T(i,j,k,l);                
               end
           end
       end
  end   
end
