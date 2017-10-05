function validation_rigid()
% validation_rigid.m
%
% validation the package for rigid ellipsoids 
% compare with Fletcher(2009) 2d
%--------------------------------------------------------------------------
   
%  clear all variables, Comment Window and figures  
   clear;
   clc;
   clf;
   
% Input parameters
   %  the bulk flow field
   L     = [0.5 1 0; 0 -0.5 0; 0 0 0];
   D     = 0.5 * (L + L');
   %  shape of the inclusion (three semi-axes of the ellipsoid)    
   a     = [10; 1; 100];
   R     = a(1)/a(2);
   
   %  orientation of the inclusion with respect to the anisotrpy
   phi   = 0:180;  % 1*n matrix
   [~,n] = size(phi);
   ang2  = [90;90];
   ang2m = repmat(ang2,1,n);
   ang   = [phi;ang2m];
   %  convert three angles from degree to radian.   
   ang_r = deg2rad(ang); 
   %  obtain the transformation matrix Q from three spherical angles, Eqns(8) 
   %  -(12) in Jiang(2007a)       
   q_n = Qvec(ang_r);
   
   m     = 25;   %  anistropy for matrix, eta_n/eta_s
   
%  generate 4th-order identity tensors   
   [Jd, ~, Ja, ~] = FourIdentity();
   
   %  viscosity of the matrix
   Cm = 2*Jd;
   Cm(1,2,:,:) = Cm(1,2,:,:)/m;
   Cm(2,1,:,:) = Cm(2,1,:,:)/m;
   Cm(2,3,:,:) = Cm(2,3,:,:)/m;
   Cm(3,2,:,:) = Cm(3,2,:,:)/m;
   
%  obtain weights and nodes before the loop 
   gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
%  allocate   
   w12D = zeros(1,n);
   for i = 1:n
      q = q_n(:,:,i); 
%  express bulk D and W in local coordinate system
      D_bar   = q * D * q';
      Cmc     = Transform(Cm,q); % Cm in clast coordinate system
%  rewrite the matrix stiffness tensor into a 1D array format 
       Carray = C2OneDarray(Cmc);
%  compute the 4th-order Green tensor T
       T      = TGreen(a, Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww);
      
%  calculate Eshelby tensors(S, PI) based on T, Eqs(3) in Qu et al.(in review)
      z       = Contract(T,Cmc);
      S       = Contract(Jd,z);
      PI      = Contract(Ja,z);
      
%  vorticity of the clast - bulk vorticity      
     invS     = FourTensorInv(S);
      u1      = Contract(PI, invS);
      wd      = Multiply(u1, D_bar);
      
      w12D(i) = -wd(1,2);
          
   end
   
% 2d Fletcher(2009)
  dxx     = D(1,1);
  dxy     = D(1,2);
  tempw12 = (R^2-1)/(R^2+1);
  phi     = deg2rad(phi);
  w12d    = -tempw12*E12(dxx,dxy,phi);
  
  phi = rad2deg(phi);
% plot
  plot(phi,w12D,'g*',phi,w12d,'b--','LineWidth',1)
  axis square
  %xlabel('\phi')
 % ylabel('w_{12}')
end


function res=E12(pure,simple,angle)
   res=-pure*sin(2*angle)+simple*cos(2*angle);
end

