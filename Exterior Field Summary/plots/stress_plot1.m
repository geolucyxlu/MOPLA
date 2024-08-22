load('InNExfields_elastic_sphere.mat', 'output');

% yz plane
stress1 = output.stress{1};
n = length(stress1);
dstress = zeros(size(stress1));
SI_d    = zeros(1,n);
I2_d    = SI_d;
I2      = SI_d;
for i=1:n
    % find deviatoric stress
    s = stress1(:,:,i);
    I2(i) = 1/2*((trace(s))^2-trace(s^2));
    m = trace(s)/3;
    s_m = m*eye(3);
    s_d = s-s_m;
    dstress(:,:,i)=s_d;
    I2_d(i) = 1/2*((trace(s_d))^2-trace(s_d^2));
    SI_d(i) = Inva(s_d);
end

pressure1 = output.pressure{1};
vonMises1 = output.vonMises{1};

X1 = squeeze(output.grid{1,1});
Y1 = squeeze(output.grid{1,2});
Z1 = squeeze(output.grid{1,3});

pressure1 = reshape(pressure1/1e9,size(Y1));
vonMises1 = reshape(vonMises1/1e9,size(Y1));
smagnitude = reshape(sqrt(I2)/1e9,size(Y1));
sdmagnitude = reshape(sqrt(-I2_d)/1e9,size(Y1));
SI_d1 = reshape(SI_d/1e9,size(Y1));

f1=figure;
subplot(2,2,1) 
contourf(Y1,Z1,pressure1)
axis equal
colormap parula
c = colorbar;
c.Label.String = 'pressure (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('pressure yz-plane')


subplot(2,2,2)
contourf(Y1,Z1,vonMises1)
axis equal
colormap parula
c = colorbar;
c.Label.String = 'von Mises stress (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('von Mises stress yz-plane')

subplot(2,2,3)
contourf(Y1,Z1,smagnitude)
axis equal
colormap parula
c = colorbar;
c.Label.String = 'stress magnitude (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('stress magnitude (I2) yz-plane')

subplot(2,2,4)
contourf(Y1,Z1,sdmagnitude)
axis equal
colormap parula
c = colorbar;
c.Label.String = 'deviatoric stress magnitude (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('deviatoric stress magnitude yz-plane')

f2=figure;
contourf(Y1,Z1,SI_d1)
axis equal
colormap parula
c = colorbar;
c.Label.String = 'SI (GPa)';
xlabel('y-axis')
ylabel('z-axis')
title('SI yz-plane')