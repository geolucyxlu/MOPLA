load('InNExfields_elastic_sphere.mat', 'output');

% yz plane
stress1 = output.stress{1};

s11 = stress1(1,1,:);
s22 = stress1(2,2,:);
s33 = stress1(3,3,:);
s23 = stress1(2,3,:);

X1 = squeeze(output.grid{1,1});
Y1 = squeeze(output.grid{1,2});
Z1 = squeeze(output.grid{1,3});

s11 = reshape(s11/1e9,size(Y1));
s22 = reshape(s22/1e9,size(Y1));
s33 = reshape(s33/1e9,size(Y1));
s23 = reshape(s23/1e9,size(Y1));

f1=figure;
subplot(2,2,1) 
contourf(Y1,Z1,s11)
axis equal
colormap parula
c = colorbar;
c.Label.String = '\sigma_{11} (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('\sigma_{11} yz-plane')


subplot(2,2,2)
contourf(Y1,Z1,s22)
axis equal
colormap parula
c = colorbar;
c.Label.String = '\sigma_{22} (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('\sigma_{22} yz-plane')

subplot(2,2,3)
contourf(Y1,Z1,s33)
axis equal
colormap parula
c = colorbar;
c.Label.String = '\sigma_{33} (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('\sigma_{33} yz-plane')

subplot(2,2,4)
contourf(Y1,Z1,s23)
axis equal
colormap parula
c = colorbar;
c.Label.String = '\sigma_{23} (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('\sigma_{23} yz-plane')
