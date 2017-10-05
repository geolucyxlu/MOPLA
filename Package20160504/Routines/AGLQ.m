function r = AGLQ(sub, x, a, b, c, d, p, ww, epsilon, C, k, varargin)
% AGLQ.m
% Adaptive Gauss-Legendre Quadrature
%
% Input:  sub, 4 subscripts of the 4th order tensor T, 1 by 4 matrix; {????}
%           x, 3 semi-axes of an inclusion, a 3 by 1 matrix;
%         a,b, the integrating range of theta; 
%         c,d, the integrating range of phi;
%     epsilon: the absolute tolerance;
%           C, stiffness of the matrix, a 1*81 matrix;

% Output:  r,  the result.
%--------------------------------------------------------------------------

  k = k + 1;
  if isempty(varargin)
     whole = dgauss(sub, x, a, b, c, d, p, ww, C);
  else
     whole = cell2mat(varargin);
  end
  n     = 0.5 * [(a+b); (c+d)];
  q(1)  = dgauss(sub, x, a, n(1), c, n(2), p, ww, C);
  q(2)  = dgauss(sub, x, a, n(1), n(2), d, p, ww, C);
  q(3)  = dgauss(sub, x, n(1), b, c, n(2), p, ww, C);
  q(4)  = dgauss(sub, x, n(1), b, n(2), d, p, ww, C);
  
  ANS    = sum(q);
  errbnd = abs(ANS - whole);
  
  if (errbnd <= max(epsilon, 10^(-3)*abs(ANS)))||(k > 100)                             
      r = whole; 
  else
      r = AGLQ(sub, x, a, n(1), c, n(2), p, ww, 0.25*epsilon, C, k, q(1)) + ...
          AGLQ(sub, x, a, n(1), n(2), d, p, ww, 0.25*epsilon, C, k, q(2)) + ...
          AGLQ(sub, x, n(1), b, c, n(2), p, ww, 0.25*epsilon, C, k, q(3)) + ...
          AGLQ(sub, x, n(1), b, n(2), d, p, ww, 0.25*epsilon, C, k, q(4));
  end    

end


function rr = dgauss(sub, x, a, b, c, d, p, ww, C)

  alpha  = 0.5*(b-a)*p + 0.5*(b+a);
  beta   = 0.5*(d-c)*p + 0.5*(d+c);
  [n,~]  = size(alpha);
  Alp    = zeros(n,n);
  Bet    = zeros(n,n);
  
  for i = 1:n
    for j = 1:n
        Alp(i,j) = alpha(i);                                          
        Bet(i,j) = beta(j);
     end
  end
  
  ww = reshape(ww,1,[]);
  x  = x';
  Alp = reshape(Alp,1,[]);
  Bet = reshape(Bet,1,[]);
  k   = partTGL(sub(1),sub(2),sub(3),sub(4),ww,Alp,Bet,C,x);
  rr  = 0.25*(d-c)*(b-a)*k;
  rr  = x(1)*x(2)*x(3)/(4*pi)*rr;
  
end

