function [ t, p ] = xyz2thetaphi ( x, y, z )

%*****************************************************************************80
%
%% XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates on the unit sphere.

%  change the range of theta from -180~180 to 0~360, 
%  and vectorize the original code, and theta, phi in rad ;
%  theta: 1*n, phi: 1*n (Mengmeng 03/24/2015)

%  Modified:
%
%    09 September 2010
%
%  Author:
%
%    Dmitri Laikov
%
%  Reference:
%
%    Vyacheslav Lebedev, Dmitri Laikov,
%    A quadrature formula for the sphere of the 131st
%    algebraic order of accuracy,
%    Russian Academy of Sciences Doklady Mathematics,
%    Volume 59, Number 3, 1999, pages 477-481.
%
%  Parameters:
%
%    Input, real X, Y, Z, the Cartesian coordinates of a point
%    on the unit sphere.
%
%    Output, real T, P, the Theta and Phi coordinates of
%    the point.
%
  p = acos ( z );
  p = p';

  fact = sqrt ( x.^2 + y.^2 );
  nonz = find(fact); % return the indexes of nonzero elements in fact.
  ze = find(fact==0); % return the indexes of zero elements in fact.
  
  t(nonz) = acos(x(nonz)./fact(nonz));
  t(ze) = acos(x(ze));
  
  id = find(y < 0.0);
  t(id) = t(id) + pi;
  
  return
end
