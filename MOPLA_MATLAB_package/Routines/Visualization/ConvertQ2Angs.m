function [ang1, ang2, ang3] = ConvertQ2Angs(q)
% ConvertQ2Angs.m
% Convert tranformation matrix Q to Theta[0,2pi] and Phi[0,pi] (i.e., cover 
% the whole range of 0<=Phi<=pi) of a1, a2, a3.
%
% Input:  q,    3*3*n matrix;

% Output: ang1,   theta, phi of a1, in radian, 2*n matrix;
%         ang2,   theta, phi of a2, in radian, 2*n matrix;
%         ang3,   theta, phi of a3, in radian, 2*n matrix;
%update:2020_4_26:converting the vectors pointing z-negtive direction 
%                 to z-positive direction 
%--------------------------------------------------------------------------

  siz = size(q);
  
  a1  = reshape(q(1,:,:),3,siz(3));
  a2  = reshape(q(2,:,:),3,siz(3));
  a3  = reshape(q(3,:,:),3,siz(3));
  
  ang1 = Vec2Sph_pos(a1);
  ang2 = Vec2Sph_pos(a2);
  ang3 = Vec2Sph_pos(a3);
end

function angs = Vec2Sph_pos(u)
% (Lecture Notes, Jiang, 2016)

% Input:  u, 3*n matrix, each col stands for a orientation vector

% Output: angs, 2*n matrix, change to theta and phi

  siz          = size(u);
  angs         = zeros(2,siz(2));
  % find the vectors pointing z-negtive direction and make them opposite 
  [~, neg]     = find(u(3,:)<0);
  u(:,neg)     = -u(:,neg);
  
  [~, col2]    = find(u(1,:));
  angs(1,col2) = atan(u(2,col2)./u(1,col2));
  [~, col3]    = find(u(1,:)==0);
  angs(1,col3) = 0.5*pi*sign(u(2,col3));
  % [~, col4]    = find(u(1,:)>=0 & u(2,:)<0);
 % angs(1,col4) = angs(1,col4) + 2*pi;
  [~, col4]    = find(u(1,:)>=0);
  angs(1,col4) = angs(1,col4);
  [~, col5]    = find(u(1,:)<0);
  angs(1,col5) = angs(1,col5) + pi;
  
  normU        = sqrt(u(1,:).^2 + u(2,:).^2 + u(3,:).^2);
  angs(2,:)    = acos(u(3,:)./normU);
end