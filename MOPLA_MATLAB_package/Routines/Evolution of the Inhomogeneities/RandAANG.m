function [a, ang]=RandAANG(a1,n)
% RandAANG.m
%
% Generate uniformly distributed random shapes of clasts as the input of
% the Mulit_rigid and Multi_deformable functions
%
% Input:  a1,  the maximun semi-axis of clasts to generate random shapes;
%         n,   number of clasts considered;

% Output: a,   random shapes: a(1,:)from 1 to a1, a(2,:)from 1 to a(1,:),
%              a(3,:)=1, 3*n matrix;
%         ang, uniform distribution of clasts, 3*n matrix;
%--------------------------------------------------------------------------
   % random shapes
   a(1,:) = 1 + (a1-1)*rand(1,n);
   a(2,:) = 1 + (a(1,:)-1).*rand(1,n);
   a(3,:) = ones(1,n);
   
   % generate a uniform distribution
   ang(1,:) = 2*pi*rand(1,n);
   ang(2,:) = acos(-1 + 2*rand(1,n));
   ang(3,:) = atan(cos(ang(2,:)).*tan(pi*rand(1,n))) + ang(1,:) + 0.5*pi;

end