function Single_deformable()
% Single_deformable.m
%
% Modeling the motion and shape evolution of a single deformable ellipsoid 
% embeded in an infinite matrix.
%
% Note:
% Default input here considering an isotropic ellipsoid embedded in a planar
% anisotropic material, however, the package is applicable for general
% anisotropic incompressible viscous materials
%--------------------------------------------------------------------------

%  clear all variables, Comment Window and figures  
   clear;
   clc;
   clf;
  
%  Input parameters:   
   L     = [0 1 0; 0 0 0; 0 0 0];  %  the bulk flow field 
   a     = [5; 3; 1];	%  3 semi-axes of the ellipsoid   
   ang   = [270; 60; 13];%  3 spherical angles (Jiang2007a, in degree)
   Nm    = 1;   %  stress exponent of the matrix
   Ne    = 3;   %  stress exponent of the ellipsoid 
   r     = 2;   %  viscosity ratio at matrix strain rate, eta_clast/eta_n
   m     = 10;   %  anistropy for matrix, eta_n/eta_s    
   epsII = 0.5; %  strain rate invariant at which ellipsoid viscosity is defiend 
   tincr = 0.01;
   steps = 1000;
   mm    = 20;    %  output steps for plotting
  
%  convert three angles from degree to radian.   
   ang_r = deg2rad(ang); 
%  decompose the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W, Eqn(3) in Jiang(2007a)     
   D     = 0.5 * (L + L');
   W     = 0.5 * (L - L');
%  obtain the transformation matrix Q from three spherical angles, Eqns(8) 
%  -(12) in Jiang(2007a)       
   q     = Q(ang_r);
  
%  generate 4th-order identity tensors   
   [Jd, ~, Ja, ~] = FourIdentity();   
%  viscosity of the matrix, Eq(12) in Qu et al.(in review)
   Cm = 2*Jd;
   Cm(1,2,:,:) = Cm(1,2,:,:)/m;
   Cm(2,1,:,:) = Cm(2,1,:,:)/m;
   Cm(2,3,:,:) = Cm(2,3,:,:)/m;
   Cm(3,2,:,:) = Cm(3,2,:,:)/m;   
%  viscosity of the clast, assuming it's isotropic here
   C_clst = 2*r*Jd;
   
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
   
%  allocate Q_evl and A_evl before the loop   
   Q_evl = zeros(3,3,steps);
   A_evl = zeros(3,steps);
   
  for k = 1:steps
%  describe D,W,C in the clast's coordinate system 
       D_bar  = q * D * q';
       W_bar  = q * W * q';
       Cmc    = Transform(Cm,q); 
       C_clstc = Transform(C_clst,q);
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
%  vorticity of the clast      
      invS    = FourTensorInv(S);
      u1      = Contract(PI, invS);
      u2      = de - D_bar;
      we      = Multiply(u1, u2) + W_bar;
      wE      = Wd(a, W_bar, de);
%  the angular velocity tensor of the ellipsoid's axes      
      Ang_vel      = we - wE;
%  update Q      
      qq           = (RodrgRot(-Ang_vel * tincr)) * q;
%  record previous Q   
      Q_evl(:,:,k) = q;
%  update a
      aa           = a.*exp(diag(de)*tincr);
%  write updated a to A_evl      
      A_evl(:,k)   = aa;
      
      qa = [qq,aa];
      qa = sortrows(qa,-4);

%     make sure that the clast is reasonable to be considered as an
%     elliposid, which means a1:a3<=100 or a2:a3<=100, boudin
      if qa(1,4)/qa(3,4) > 100
          qa(1,4) = 0.5*qa(1,4);
      elseif qa(2,4)/qa(3,4) > 100
          qa(2,4) = 0.5*qa(2,4);
      else
      end 
      qa = sortrows(qa,-4);
      a  = qa(:,4);       
      q  = qa(1:3,1:3);
  end
  
% save Q_evl and A_evl to the current workspace 
 % save('evl_single_deformable.mat','Q_evl','A_evl');
  
% Output steps for plotting the rotation path
   nn       = 1:mm:steps; 
   [~,last] = size(nn);
   Q_plot   = Q_evl(:,:,nn);
   A_plot   = A_evl(:,nn);
   
   
%  Equal-area projection
   
%     compute two spherical angles for three axes
      [a1_evl, a2_evl, a3_evl] = ConvertQ2Angs(Q_plot);
   
%     compute r for equal-area projection, both hemispheres will be plotted
%     a1
      [~,a1in]  = find(a1_evl(2,:)<=(0.5*pi));
      [~,a1out] = find(a1_evl(2,:)>(0.5*pi));
      r1(a1in)  = sqrt(2) * sin(a1_evl(2,a1in)./2);
      r1(a1out) = sqrt(2) * cos(a1_evl(2,a1out)./2);
%     a2
      [~,a2in]  = find(a2_evl(2,:)<=(0.5*pi));
      [~,a2out] = find(a2_evl(2,:)>(0.5*pi));
      r2(a2in)  = sqrt(2) * sin(a2_evl(2,a2in)./2);
      r2(a2out) = sqrt(2) * cos(a2_evl(2,a2out)./2);
%     a3       
      [~,a3in]  = find(a3_evl(2,:)<=(0.5*pi));
      [~,a3out] = find(a3_evl(2,:)>(0.5*pi));
      r3(a3in)  = sqrt(2) * sin(a3_evl(2,a3in)./2);
      r3(a3out) = sqrt(2) * cos(a3_evl(2,a3out)./2);
      
%     equal-area projections of a1, a2, a3   
%     a1
   
      subplot(1,3,1);
      t = 0 : .01 : 2 * pi;
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
%     phi<=pi/2, plot red dots 
      polar(a1_evl(1,a1in),r1(a1in),'.r')
%     phi>pi/2, plot green dots 
      polar(a1_evl(1,a1out),r1(a1out),'.g')
%     starting point      
      polar(a1_evl(1,1),r1(1),'xb')
%     end point 
      polar(a1_evl(1,last),r1(last),'*c')
      hold off
      title('a1')
      
%     a2   
      subplot(1,3,2);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a2_evl(1,a2in),r2(a2in),'.r')
      polar(a2_evl(1,a2out),r2(a2out),'.g')
      polar(a2_evl(1,1),r2(1),'xb')
      polar(a2_evl(1,last),r2(last),'*c')
      hold off
      title('a2')
      
%     a3   
      subplot(1,3,3);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a3_evl(1,a3in),r3(a3in),'.r')
      polar(a3_evl(1,a3out),r3(a3out),'.g')
      polar(a3_evl(1,1),r3(1),'xb')
      polar(a3_evl(1,last),r3(last),'*c')
      hold off
      title('a3')
   
   
  
   % Flinn diagram
    x = log(A_plot(2,:)./A_plot(3,:));
    y = log(A_plot(1,:)./A_plot(2,:));
    
    figure('Name','Flinn diagram: The shape evolution of a deformable ellipsoid');
    plot(x,y,'.r',0:3,0:3,'-k')
    hold on
    plot(x(1),y(1),'xc','MarkerSize',5)
    plot(x(last),y(last),'*b','MarkerSize',5)
    xlabel('ln(a2/a3)')
    ylabel('ln(a1/a2)')
    axis square
  
end