function [de,C_clst] = Ed(Nm,Nc,Cm,Ce,r,s,d,epsilonII,Jd,q)
% Ed.m
%
% Iterative calculation in the ellipsoid, starting with current strain rate
% invariant "REF", until the new strain rate invariant converges. 
%
% Input: Nm,        power-law stress exponent of the matrix;
%        Nc,        power-law stress exponent of the ellipsoid;
%        r,         the viscosity ratio(at matrix strain rate); 
%        s,         the Eshelby tensor S, a 4th order tensor;
%        d,         the imposed strain rate tensor expressed in x' system 
%                   at current state, 3*3 matrix;
%        epsilonII, the strain rate invariant at which ellipsoid viscosity 
%                   is defined;

% Output: de,       the strain rate tensor of the ellipsoid, 3*3 matrix;
%         zz,       the updated r.
%--------------------------------------------------------------------------
  yy      = r;
  REF     = epsilonII;
  C_clst  = 2*yy*Jd; 
  
  for jj = 1:1000
      f        = fdE(Nm,Cm,Ce,s,Jd);
      de       = Multiply(f,d);
%     the second invariant of the strain rate tensor
      epsilonI = (0.5*contract1(de,de))^0.5;
      zz       = (epsilonI/REF)^((1-Nc)/Nc)*yy;
      alpha    = abs(zz/yy - 1);
      yy       = zz;
      REF      = epsilonI;
      if alpha < 0.001
          break
      end
	  C_clst  = 2*yy*Jd; 
      Ce = Transform(C_clst,q);
  end
  
end

function r = fdE(n,Cm,Ce,S,Jd)
      Mm = FourTensorInv(Cm);
      bb = Contract(S,Contract(Mm,Ce));
      b = Jd-S+n*bb;
      a = FourTensorInv(b);
      c = Jd + (n-1)*S;                                                  
      r = Contract(a,c);
end
