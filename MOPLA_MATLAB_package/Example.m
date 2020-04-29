% Cascade Lake shear zone (Jiang, 2014)

% clear all variables, Comment Window and figures--------------------------  
   clear;
   clc;
   
%% Imposed macroscale flow filed
   alpha = deg2rad(15);
   beta  = deg2rad(70);
   L     = [0  -cos(alpha)          0;...
            0 -sin(alpha)*sin(beta) 0;...
            0  sin(alpha)*cos(beta) sin(alpha)*sin(beta)];
%% Primary inclusions set up 
% The number of primary inclusions.
  n = 200;
% Generate uniformly distributed primary inclusions with random shapes.
  [a, ang] = RandAANG(10,n);
% Generate the transformation tensors(q) of primary inclusions.
  q  = Qvec(ang);  
% Generate the effective viscosities(eta) of primary inclusions 
% at the macroscale strain-rate state (the reference state).
  eta(1,1:120) = ones(1,120);
  vmin = 1;
  vmax = 10;
  eta(1,121:200) = vmin + (vmax-vmin)*rand(1,80);
% Generate the stress exponents(Ne) of primary inclusions.
  Nmax = 4;
  Nmin = 2;
  Ne   = Nmin + (Nmax-Nmin)*rand(1,n);
%% Secondary inclusions set up 
% The number of secondary inclusion in each primary inclusion.
  ns   = 2;
% 1st secondary-order insluion 
% (uniformly distributed, with random shapes)
 [ak(:,:,1),angk1] = RandAANG(10,n);
 qk(:,:,:,1) = Qvec(angk1);
 etak(1,:)= 2.*eta;
 Nk(1,:)= 3.* ones(1,n);
% 2nd secondary-order insluion 
 ak(:,:,2) = ones(3,n);
 [~,angk2] = RandAANG(10,n);
 qk(:,:,:,2) = Qvec(angk2);
 etak(2,:) = eta;
 Nk(2,:) = Ne;
%% Modeling the mechanicall behavior of the heterogeneous rock mass 
% The time increment of each computational step.  
  tincr = 0.01;
% total computational steps.
  steps = 250;
   
% % MOPLA Model with primary inclusions  
%  [S_bar_evl,C_bar_evl,Q_evl,A_evl]= MOPLA_primary(L, n, a, q, eta, Ne, steps, tincr);
[Sm_evl,S_evl,C_bar_evl,Q_evl,A_evl,Qk_evl,Ak_evl]= MOPLA_secondary(L, n, a, q, eta, Ne, steps, tincr, ns, ak, qk, etak, Nk);

%% Visualization 
% %  Plots of primary ellipsoids
% %  Get the final orientaions of the primary ellipsoids
%    Q_final = Q_evl(:,:,:,steps); 
%    f1=figure('Name','Stereonet: The orientation of ellipsoids');
%    Stereonet(Q_final);
%    
% %  Get the final shape of the primary ellipsoids   
%    A = A_evl(:,:,steps);
%    f2=figure('Name','Flinn diagram: The shape of ellipsoids');
%    Flinn(A);
   
  %  Plots of secondary ellipsoids 
   Qk_1 = Qk_evl(:,:,:,1,steps); 
   f3=figure('Name','Stereonet: The orientation of NO.1 secondary inclusion');
   Stereonet(Qk_1);
   
   Ak_1 = Ak_evl(:,:,1,steps); 
   f4=figure('Name','Flinn diagram: The shape of NO.1 secondary inclusion');
   Flinn(Ak_1);  
   
   
   Qk_2 = Qk_evl(:,:,:,2,steps); 
   f5=figure('Name','Stereonet: The orientation of NO.2 secondary inclusion');
   Stereonet(Qk_2);
   
   Ak_2 = Ak_evl(:,:,2,steps); 
   f6=figure('Name','Flinn diagram: The shape of NO.2 secondary inclusion');
   Flinn(Ak_2);  
   
