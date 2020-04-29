function [p,w] = Gauss(n)
% Gauss.m
% Gauss-Legendre nodes and weights
%
% Input:  n, the number of points wanted to generate;

% Output: p, nodes, n*1 matrix;
%         w, weights, n*1 matrix.
%--------------------------------------------------------------------------
  x = zeros(n,1);
  w = zeros(n,1);
  m = 0.5*(n+1);
  
  for ii = 1:m
     z  = cos(pi*(ii-0.25)/(n+0.5));
     z1 = z + 1;
     while abs(z-z1) > 10^(-10)
        p1 = 1;
        p2 = 0;
        for jj = 1:n
           p3 = p2;
           p2 = p1;
           p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj;
        end
        pp = n*(z*p1-p2)/(z^2-1);
        z1 = z;
        z  = z1 - p1/pp;
     end
     x(ii,1)     = -z;
     x(n+1-ii,1) = z;
     w(ii,1)     = 2/((1-z^2)*(pp^2));
     w(n+1-ii,1) = w(ii,1);
  end
  
  p = x;
end