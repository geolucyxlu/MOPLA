%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the data structure needed for an ellipsoid inhomogeniety problem;
% solve the problem with an equevilent Eshelby inclusion problem;
% demonstrate the outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Sharma's ellpdoidal heterogeneity with observation line
%clear all;
 K0 = 156e9;
 G0 = 209e9;
 K1 = 83e9;  
 G1 = 107e9;
 nu0 = (3*K0 - 2*G0)/2/(3*K0 + G0);
 nu1 = (3*K1 - 2*G1)/2/(3*K1 + G1);
 E0  = 9*K0*G0/(3*K0+G0);
 E1  = 9*K1*G1/(3*K1+G1);
%poisson ratio of matrix
%incl.vm= 0.3465;
incl.vm= nu0;
%elastic modulus of matrix
%incl.Em=70e9;
incl.Em=E0;
%hetergeneity poisson ratio
%incl.vh=0.2324;
incl.vh=nu1;
%hetergeneity elastic modulus
%incl.Eh=106e9;
incl.Eh=E1;

% dimensiona of the ellipsoid.
incl.dim=[10000 1.000001 1];
% ortation angles in radian applied to the ellipsoid in sequence of [x y z];.
incl.ang = [0 0 0];

% remote stress ordered by sigma11, sigma12, sigma13, sigma22,sigma23,sigma33
incl.stressvec=[-5e9;0;0;-8e9;0;-2e9];
incl.eigp=[-0.05;0;0;-0.05;0;-0.05];
% first observaton grid
x = 0;
y = -5:.01:5;
z = -5:.01:5;
incl.grid{1} = {x,y,z};
% % second observation grid
% x = -5:.01:5;
% y = 0;
% z = -5:.01:5;
% incl.grid{2} = {x,y,z};
% % third observation grid
% x = -5:.01:5;
% y = -5:.01:5;
% z = 0;
% incl.grid{3} = {x,y,z};
clear x y z
% add more grids if needed

% call Eshelby solver,arg: 'disp','stress','strain' only output displacements by default
incl.sol = Esh_sol(incl,'disp','stress','strain');

% demonstrate the third grid's results for the first and second components of
% the output displacement and stress.
X = squeeze(incl.sol.grid{1,1});
Y = squeeze(incl.sol.grid{1,2});
Z = squeeze(incl.sol.grid{1,3});
u = squeeze(incl.sol.u{1});
stress = squeeze(incl.sol.stress{1});

save('normalstress3GPa_1.mat', 'incl')
% az = 90;
% el = 0;
% subplot(2,2,1)
% surf(X,Y,stress(:,:,1),'LineStyle','none');
% colorbar;
% view(az, el);
% subplot(2,2,2)
% surf(X,Y,stress(:,:,4),'LineStyle','none');
% colorbar;
% view(az, el);
% subplot(2,2,3)
% surf(X,Y,stress(:,:,6),'LineStyle','none');
% colorbar;
% view(az, el);
% subplot(2,2,4)
% surf(X,Y,stress(:,:,4),'LineStyle','none');
% view(az, el);
% colorbar;

