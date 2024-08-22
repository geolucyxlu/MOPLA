 %% plot
 load('InNExfields_elastic_sphere.mat','input','output')
  % first observation grid
  pressure1 = output.pressure{1};
  vonMises1 = output.vonMises{1};

  X1 = squeeze(output.grid{1,1});
  Y1 = squeeze(output.grid{1,2});
  Z1 = squeeze(output.grid{1,3});

  pressure1 = reshape(pressure1/1e9,size(Y1));
  vonMises1 = reshape(vonMises1/1e9,size(Y1));

  % second observation grid
  pressure2 = output.pressure{2};
  vonMises2 = output.vonMises{2};

  X2 = squeeze(output.grid{2,1});
  Y2 = squeeze(output.grid{2,2});
  Z2 = squeeze(output.grid{2,3});

  pressure2 = reshape(pressure2/1e9,size(X2));
  vonMises2 = reshape(vonMises2/1e9,size(X2));

  % third observation grid
  pressure3 = output.pressure{3};
  vonMises3 = output.vonMises{3};

  X3 = squeeze(output.grid{3,1});
  Y3 = squeeze(output.grid{3,2});
  Z3 = squeeze(output.grid{3,3});

  pressure3 = reshape(pressure3/1e9,size(X3));
  vonMises3 = reshape(vonMises3/1e9,size(X3));
 

  % Pressure plot
  % figure('Name','Stress plot');
   subplot(3,2,1) 
   contourf(Y1,Z1,pressure1)
   axis equal
   colormap parula
   c = colorbar;
   c.Label.String = 'pressure (GPa)';
   xlabel('y-axis')
   ylabel('z-axis')
   title('pressure (yz plane)')

   subplot(3,2,2)
   contourf(Y1,Z1,vonMises1)
   axis equal
   colormap parula
   c = colorbar;
   c.Label.String = 'von Mises stress (GPa)';
   xlabel('y-axis')
   ylabel('z-axis')
   title('von Mises stress (yz plane)')

   subplot(3,2,3)
   contourf(X2,Y2,pressure2)
   axis equal
   colormap parula
   c = colorbar;
   c.Label.String = 'pressure (GPa)';
   xlabel('x-axis')
   ylabel('y-axis')
   title('pressure (xy plane)')

   subplot(3,2,4)
   contourf(X2,Y2,vonMises2)
   axis equal
   colormap parula
   c = colorbar;
   c.Label.String = 'von Mises stress (GPa)';
   xlabel('x-axis')
   ylabel('y-axis')
   title('von Mises stress (xy plane)')

   subplot(3,2,5)
   contourf(X3,Z3,pressure3)
   axis equal
   colormap parula
   c = colorbar;
   c.Label.String = 'pressure (GPa)';
   xlabel('x-axis')
   ylabel('z-axis')
   title('pressure (xz plane)')

   subplot(3,2,6)
   contourf(X3,Z3,vonMises3)
   axis equal
   colormap parula
   c = colorbar;
   c.Label.String = 'von Mises stress (GPa)';
   xlabel('x-axis')
   ylabel('z-axis')
   title('von Mises stress (xz plane)')