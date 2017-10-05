function [PHI, R] = Incl(a)
% Incl.m
% get PHI, R of an ellipsoid
%
% PHI measures of the triaxiality of an ellipsoid
% R measures of ellipsoid deviation from sphere

% Input:  a are three semi-axes of n ellipsoids, 3*n matrix

% Output: PHI, 1*n matrix, in deg
%         R,   1*n matrix
%--------------------------------------------------------------------------
  
  a = sort(a,'descend');
  x = log(a(2,:)./a(3,:));
  y = log(a(1,:)./a(2,:));
  
  nor       = find(x~=0); 
  k         = y(nor)./x(nor);  
  PHI(x==0) = 90; 
  PHI(nor)  = atan(k)*180/pi;
 
  R   = sqrt(y.^2 + x.^2);

end