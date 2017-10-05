function Multi_deformable()
% Multi_deformable.m
%
% Modeling the evolution of a population of n dilute (non-interacting) defomrable  
% ellipsoids embedded in a general viscous material undergoing homogeneous deformation.
%
% Note:
% Default input here considering isotropic ellipsoids embedded in a planar
% anisotropic material, however, the package is applicable for general
% anisotropic incompressible viscous materials
%--------------------------------------------------------------------------

%  clear all variables, Comment Window and figures  
   clear;
   clc;
   clf;
   
% Input parameters
% the bulk flow field
   L     = [1 0 0; 0 -1 0; 0 0 0];
%  number of clasts considered;   
   n     = 100;
%  the maximun semi-axis of clasts to generate random shapes; 
   a1    = 10;
%  stress exponent of the matrix  
   Nm    = 1;
%  stress exponent of the ellipsoid
   Ne    = 3;
%  maximum viscosity ratio at matrix strain rate    
   rmax  = 3;
%  minimum viscosity ratio at matrix strain rate       
   rmin  = 1;
%  strain rate invariant at which ellipsoid viscosity is defiend  
   epsII = 0.5;
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.01;
%  total steps of the computation
   steps = 500;
%  anistropy for matrix, eta_n/eta_s          
   m     = 5;   
%  Initial shapes and orientations of clasts:
%  generate a population of uniformly distributed ellipsoids with
%  a1:a2:1,a2 from a1 to 1.
   [a, ang] = RandAANG(a1,n);
  
% input or generation of a population of initial effective viscosity ratios r
% eta_clast/eta_n
   r = rmin + (rmax-rmin)*rand(1,n);
 
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W, Eqn(3) in Jiang(2007a)        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');
   
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

   % allocate variables   
   Q_evl = zeros(3,3,steps,n);
   A_evl = zeros(3,steps,n);
  
   for l = 1:n 
       q      = Q(ang(:,l));
       C_clst = 2*r(l)*Jd;
       for k = 1:steps
%  describe D,W,C in the clast's coordinate system 
       D_bar   = q * D * q';
       W_bar   = q * W * q';
       Cmc     = Transform(Cm,q); 
       C_clstc = Transform(C_clst,q);
%  rewrite the matrix stiffness tensor into a 1D array format 
       Carray  = C2OneDarray(Cmc);
%  compute the 4th-order Green tensor T
       T       = TGreen(a(:,l), Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww);  
           
%  calculate Eshelby tensors(S, PI) based on T, Eqs(3) in Qu et al.(in review)
      z        = Contract(T,Cmc);
      S        = Contract(Jd,z);
      PI       = Contract(Ja,z);
%  strain rate of the clast      
      [de, C_clst] = Ed(Nm, Ne, Cmc, C_clstc, r(l), S, D_bar, epsII, Jd, q);    
%  vorticity of the clast      
      invS    = FourTensorInv(S);
      u1      = Contract(PI, invS);
      u2      = de - D_bar;
      we      = Multiply(u1, u2) + W_bar;
      wE      = Wd(a(:,l), W_bar, de);
%  the angular velocity tensor of the ellipsoid's axes      
      Ang_vel      = we - wE;
%  update Q      
      qq           = (RodrgRot(-Ang_vel * tincr)) * q;
%  record previous Q   
		 Q_evl(:,:,k,l) = q;
%     update a         
         aa             = a(:,l).*exp(diag(de)*tincr);
%     write updated a to A_evl         
         A_evl(:,k,l)   = a(:,l);
%     make sure that Q and a are in the descending oreder of a(a1>=a2>=a3) 
         qa = [qq,aa];
         qa = sortrows(qa,-4);
         
%     ensure that the clast is not too elongated or flattened 
%     setting a1:a3<=100 or a2:a3<=100, boudinage otherwise
          if qa(1,4)/qa(3,4) > 100
              qa(1,4) = 0.5*qa(1,4);
          elseif qa(2,4)/qa(3,4) > 100
              qa(2,4) = 0.5*qa(2,4);
          else
          end          
         qa         = sortrows(qa,-4);     
         a(:,l)     = qa(:,4);       
         q = qa(1:3,1:3);
       end
   end
   
% save Q_evl and A_evl to the current workspace   
  %save('evl_multi_deformable.mat','Q_evl','A_evl');
   
%  get the final orientaions of the rigid ellipsoids
   Q_final = squeeze(Q_evl(:,:,steps,:));

%  Equal-area projection

%     compute two spherical angles for three axes
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
   
% Flinn diagram
    A = squeeze(A_evl(:,steps,:));
    x = log(A(2,:)./A(3,:));
    y = log(A(1,:)./A(2,:));
    
    figure('Name','Flinn diagram: The shape of deformable ellipsoids');
    
    plot(x,y,'.r',0:5,0:5,'-k')
    xlabel('ln(a2/a3)')
    ylabel('ln(a1/a2)')
    axis square
end
