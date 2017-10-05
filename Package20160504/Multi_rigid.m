function Multi_rigid()
% Multi_rigid.m
%
% Modelling the rotation and preferred orientation development for a system
% of multiple non-interacting rigid ellipsoids embedded in general 
% anisotropic incompressible viscous matrix
%
%--------------------------------------------------------------------------
 
%  clear all variables, Comment Window and figures  
   clear;
   clc;
   clf;
   
%  Input parameters:    
   L     = [0.375 1 0; 0 -0.375 0; 0 0 0];  %  the bulk flow field 
   n     = 300;   % the number of the ellipsoids
   a1    = 10;     % the longest semi-axis of those ellipsoids
   tincr = 0.05;  % the time increment
   steps = 200;  % total steps    
   m     = 1;    %  anistropy for matrix, eta_n/eta_s
   
%  Initial shapes and orientations of clasts:
%  generate a population of uniformly distributed ellipsoids with
%  a1:a2:1,a2 from a1 to 1.
   [a, ang] = RandAANG(a1,n);
   
%  decompose the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W, Eqn(3) in Jiang(2007a)      
   D     = 0.5 * (L + L');
   W     = 0.5 * (L - L');
%  generate 4th-order identity tensors   
   [Jd, ~, Ja, ~] = FourIdentity();
%  viscosity of the matrix, Eq(12) in Qu et al.(in review)
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
   
%  allocate Q_evl before the loop     
   Q_evl = zeros(3,3,steps,n);
   
%  start calculating the rotation of the inclusion, Eqs(9) in Qu et al.(in review)
%  applying Rodrigues' rotation approximation(Jiang, 2013) to solve Eq(9a) 
   for l = 1:n     
       q = Q(ang(:,l));        
       for k = 1:steps
%  describe D,W,C in the clast's coordinate system 
       D_bar  = q * D * q';
       W_bar  = q * W * q';
       Cmc    = Transform(Cm,q); 
%  rewrite the matrix stiffness tensor into a 1D array format 
       Carray = C2OneDarray(Cmc);
%  compute the 4th-order Green tensor T
       T      = TGreen(a(:,l), Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww);           
%  calculate Eshelby tensors(S, PI) based on T, Eqs(3) in Qu et al.(in review)
      z       = Contract(T,Cmc);
      S       = Contract(Jd,z);
      PI      = Contract(Ja,z);
%  update the angular velocity of the ellipsoid, Eq(9c) in Qu et al.(in review)       
      invS    = FourTensorInv(S);
      u1      = Contract(PI, invS);
      wd      = Multiply(u1, D_bar);
      Ang_vel = W_bar - wd;
%  Rodrigues' rotation approximation to update Q, Eq(40) in Jiang(2013)
       qq           = (RodrgRot(-Ang_vel * tincr)) * q;
%  record Q to Q_evl       
       Q_evl(:,:,k,l) = q;
       q              = qq;
       end       
   end
%  save Q_evl to the current workspace      
   %save('Q_evl_multi_rigid.mat','Q_evl');

   Q_final = squeeze(Q_evl(:,:,steps,:));

%  Equal-area projection

      [a1_ang, a2_ang, a3_ang] = ConvertQ2Angs(Q_final);
      
%     compute r for equal-area projection, both hemispheres will be plotted
%     a1
      [~,a1in]  = find(a1_ang(2,:)<=(0.5*pi));
      [~,a1out] = find(a1_ang(2,:)>(0.5*pi));
      r1(a1in)  = sqrt(2) * sin(a1_ang(2,a1in)./2);
      r1(a1out) = sqrt(2) * cos(a1_ang(2,a1out)./2);
%     a2
      [~,a2in]  = find(a2_ang(2,:)<=(0.5*pi));
      [~,a2out] = find(a2_ang(2,:)>(0.5*pi));
      r2(a2in)  = sqrt(2) * sin(a2_ang(2,a2in)./2);
      r2(a2out) = sqrt(2) * cos(a2_ang(2,a2out)./2);
%     a3       
      [~,a3in]  = find(a3_ang(2,:)<=(0.5*pi));
      [~,a3out] = find(a3_ang(2,:)>(0.5*pi));
      r3(a3in)  = sqrt(2) * sin(a3_ang(2,a3in)./2);
      r3(a3out) = sqrt(2) * cos(a3_ang(2,a3out)./2);
   
%     equal-area projections of a1, a2, a3   
%     a1      
      subplot(1,3,1);
      t = 0 : .01 : 2 * pi;
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
%     phi<=pi/2, plot red dots 
      polar(a1_ang(1,a1in),r1(a1in),'.r')
%     phi>pi/2, plot red dots 
      polar(a1_ang(1,a1out),r1(a1out),'.r')
      hold off
      title('a1')
      
%     a2   
      subplot(1,3,2);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a2_ang(1,a2in),r2(a2in),'.r')
      polar(a2_ang(1,a2out),r2(a2out),'.r')
      hold off
      title('a2')
      
%     a3   
      subplot(1,3,3);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a3_ang(1,a3in),r3(a3in),'.r')
      polar(a3_ang(1,a3out),r3(a3out),'.r')
      hold off
      title('a3')
end