function [r0,r1,r2,r3] = EmprFun(PHI)
% EmprFun.m
% Empirical function 
%
% Input:  PHI
%         
% Output: r0,r1,r2,r3
%--------------------------------------------------------------------------
  x = PHI;
  r0 = fitobj0(x);
  r1 = fitobj1(x);
  r2 = fitobj2(x);
  r3 = fitobj3(x);
 
end
function r = fitobj0(x)

       p1 =  -1.944e-08;
       p2 =  -4.289e-06;
       p3 =    0.001325;
       p4 =     -0.0958;
       p5 =       2.552; 
   r = p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x + p5;         
end

function r = fitobj1(x)
% cputime~15
       p1 =   1.025e-08;
       p2 =  -2.461e-06;
       p3 =   0.0006203;
       p4 =    -0.06954;
       p5 =       4.018;

r = p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x + p5;   
end

function r = fitobj2(x)
% cputime~30
       p1 =   1.496e-07;
       p2 =  -2.729e-05;
       p3 =    0.001792;
       p4 =    -0.07242;
       p5 =       4.583;
       
 r = p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x + p5;               
end

function r = fitobj3(x)
% cputime~40
       a0 =       3.693;
       a1 =      0.7857;
       b1 =      0.4061;
       a2 =      0.1859;
       b2 =    -0.09161;
       a3 =     0.08343;
       b3 =    -0.02982;
       w =     0.05009;
       
   r = a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + ...
       a3*cos(3*x*w) + b3*sin(3*x*w);
              
end

