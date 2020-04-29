function [Alp,Bet,ww] = GaussGGLQ(n)
% GaussGGLQ.m
% Gauss-Legendre method to generate nodes and weights for GGLQ
%
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
  
  ww    = w*w';
  alpha = 0.5*(2*pi-0)*p + 0.5*(2*pi+0);
  beta  = 0.5*(pi-0)*p + 0.5*(pi+0);
  [n,~] = size(alpha);
  Alp   = zeros(n,n);
  Bet   = zeros(n,n);

  for i = 1:n
    for j = 1:n
        Alp(i,j) = alpha(i);                                          
        Bet(i,j) = beta(j);
     end
  end
end