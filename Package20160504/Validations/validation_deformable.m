function validation_deformable()
% validation_deformable.m
%
% validation the package for deformable ellipsoids 
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
   a     = [5; 1; 100];
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
   
   r     = 10;   %  viscosity ratio at matrix strain rate, eta_clast/eta_n
   m     = 25;   %  anistropy for matrix, eta_n/eta_s
   
%  generate 4th-order identity tensors   
   [Jd, ~, Ja, ~] = FourIdentity();
   
   %  viscosity of the matrix
   Cm = 2*Jd;
   Cm(1,2,:,:) = Cm(1,2,:,:)/m;
   Cm(2,1,:,:) = Cm(2,1,:,:)/m;
   Cm(2,3,:,:) = Cm(2,3,:,:)/m;
   Cm(3,2,:,:) = Cm(3,2,:,:)/m;
   
   %  clast viscosity
   C_clst = 2*r*Jd;
   %  stress exponent of the matrix  
   Nm    = 1;
   %  stress exponent of the ellipsoid
   Ne    = 1;
   %  strain rate invariant at which ellipsoid viscosity is defiend  
   epsII = 0.5;
   
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
   d11D = zeros(1,n);
   d12D = zeros(1,n);
   w12D = zeros(1,n);
   for i = 1:n
      q = q_n(:,:,i); 
%  express bulk D and W in local coordinate system
      D_bar   = q * D * q';
      Cmc     = Transform(Cm,q); % Cm in clast coordinate system
      C_clstc = Transform(C_clst,q); % clast viscosity in clast coordinate system 
%  rewrite the matrix stiffness tensor into a 1D array format 
       Carray = C2OneDarray(Cmc);
%  compute the 4th-order Green tensor T
       T      = TGreen(a, Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww);
      
%  calculate Eshelby tensors(S, PI) based on T, Eqs(3) in Qu et al.(in review)
      z       = Contract(T,Cmc);
      S       = Contract(Jd,z);
      PI      = Contract(Ja,z);
%  strain rate of the clast      
      [de, C_clst]  = Ed(Nm, Ne, Cmc, C_clstc, r, S, D_bar, epsII, Jd, q);    
      
      deD     = de-D_bar;
      d11D(i) = deD(1,1);
      d12D(i) = deD(1,2);
      
%  vorticity of the clast - bulk vorticity      
      invS = FourTensorInv(S);
      u1   = Contract(PI, invS);
      u2   = de - D_bar;
      wD   = Multiply(u1, u2);
      
      w12D(i) = wD(1,2);
          
   end
   
% 2d Fletcher(2009)
  dxx     = D(1,1);
  dxy     = D(1,2);
  temp11  = 2*(m^0.5)*R/(R^2+2*r*R*(m^0.5)+1);
  temp12  = (m^0.5)*(R^2+1)/(r*(m^0.5)*(R^2+1)+2*R);
  tempw12 = (R^2-1)/(R^2+1);
  phi     = deg2rad(phi);
  
  e11d = temp11*(-r*E11(dxx,dxy,phi)+E11(dxx,dxy/m,phi));
  e12d = temp12*(-r*E12(dxx,dxy,phi)+E12(dxx,dxy/m,phi));
  w12d = tempw12*e12d;
  
  phi = rad2deg(phi);
  % plot
  subplot(1,3,1)
  plot(phi,d11D,'rx',phi,e11d,'b--','LineWidth',1)
  axis square
  
  subplot(1,3,2)
  plot(phi,d12D,'mx',phi,e12d,'b--','LineWidth',1) 
  axis square
  subplot(1,3,3)
  plot(phi,w12D,'gx',phi,w12d,'b--','LineWidth',1)
  axis square
end

function re=E11(pure,simple,angle)
    re=pure*cos(2*angle)+simple*sin(2*angle);
end

function res=E12(pure,simple,angle)
   res=-pure*sin(2*angle)+simple*cos(2*angle);
end

