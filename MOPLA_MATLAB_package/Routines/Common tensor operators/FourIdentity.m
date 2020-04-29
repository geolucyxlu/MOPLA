function [Jd, Js, Ja, Jm] = FourIdentity()
% FourIdentity.m
% 4th-order identity tensors
%
% Output: Jd, 4*4*4*4 matrix
%         Js, 4*4*4*4 matrix
%         Ja, 4*4*4*4 matrix
%         Jm, 4*4*4*4 matrix
%--------------------------------------------------------------------------
%  delta
   delta = [1 0 0; 0 1 0; 0 0 1];
  
%  allocate J, Jm
   J  = zeros(3,3,3,3);
   J1 = zeros(3,3,3,3);
   Jm = zeros(3,3,3,3);

%  Eqn(3) in Jiang(2014)
   for i = 1:3
       for j = 1:3
           for k = 1:3
               for l = 1:3
                   J(i,j,k,l)  = delta(i,k)*delta(j,l);
                   J1(i,j,k,l) = delta(j,k)*delta(i,l);
                   Jm(i,j,k,l) = delta(i,j)*delta(k,l)/3;
               end
           end
       end
   end
   Js = 0.5*(J + J1);
   Ja = 0.5*(J - J1);
   Jd = Js - Jm;   
end