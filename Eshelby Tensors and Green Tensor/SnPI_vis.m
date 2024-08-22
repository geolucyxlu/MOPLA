function [S,PI,LAMBDA ] = SnPI_vis(a)
% Calculation of interior Eshelby tensors for viscous isotropic materials
% from corresponding elastic expressions by setting the Poisson's ratio
% to 0.5.

    % Elastic Eshelby tensors (s,p) when Poisson's ratio equals to 0.5. 
    s = zeros(3,3,3,3);
    p = zeros(size(s));
    JJ = J(a);
    h = JJ(1:3,4);
    for lambda = 1:3
       for eta = 1:3
           s(lambda,eta,lambda,eta) = 3/(8*pi)*(a(lambda)^2+a(eta)^2)*JJ(lambda,eta);
           s(lambda,eta,eta,lambda) = s(lambda,eta,lambda,eta);
           s(lambda,lambda,eta,eta) = 3/(4*pi)*a(eta)^2*JJ(lambda,eta);

           p(lambda,eta,lambda,eta) = ( h(eta)-h(lambda) )/(8*pi);
           p(lambda,eta,eta,lambda) = p(lambda,eta,lambda,eta);
           p(eta,lambda,lambda,eta) = -p(lambda,eta,lambda,eta);
           p(eta,lambda,eta,lambda) = -p(lambda,eta,lambda,eta);
       end
    end
    
    % The modification for the diagonal components (S,PI,LAMDA) (Jiang, 2016)
    % LAMBDA
      LAMBDA   = zeros(3,3);
      for i=1:3
          LAMBDA(i,i) = -1/3* (s(i,i,1,1)+ s(i,i,2,2)+ s(i,i,3,3)); 
      end
    % S  
      S       = s;
      for i=1:3
          for j=1:3
              S(i,i,j,j) = LAMBDA(i,i)+ s(i,i,j,j);
          end
      end
    % PI  
      PI      = p;
end

function J = J(a)
    % J(X,Y,Z)
    if (a(1)>a(2))&&(a(2)>a(3))
        J = triaxial(a);
    elseif (a(1)==a(2))&&(a(2)>a(3))
        J = oblate(a);
    elseif (a(1)>a(2))&&(a(2)==a(3))
        J = prolate(a);
    elseif (a(1)==a(2))&&(a(2)==a(3))
        J = sphere(a);
    else
        print('wrong');
    end
end

function triaxial = triaxial(a)
    % Triaxial
    phi= Phi(a);
    V = a(1)*a(2)*a(3);
    bb = [1/(a(1)^2-a(3)^2)^0.5; 1/(a(2)^2-a(3)^2)^0.5; 1/(a(1)^2-a(2)^2)^0.5];
    I1 = 4*pi*V*bb(1)*bb(3)^2*(phi(1)-phi(2));
    I3 = 4*pi*V*bb(1)*bb(2)^2*(a(2)*(a(1)^2-a(3)^2)^0.5/(a(1)*a(3))-phi(2));
    I2 = 4*pi-I1-I3;
    I12 = (I2-I1)*bb(3)^2 /3;
    I13 = (I3-I1)*bb(1)^2 /3;
    I23 = (I3-I2)*bb(2)^2 /3;
    I11 = 4*pi/(3*a(1)^2)-I12-I13;
    I22 = 4*pi/(3*a(2)^2)-I12-I23;
    I33 = 4*pi/(3*a(3)^2)-I13-I23;
    triaxial = [I11 I12 I13 I1;I12 I22 I23 I2;I13 I23 I33 I3];
end

function oblate = oblate(a)
    % Oblate
    bb1 = a(1)^2-a(3)^2;
    bb2 = a(3)/a(1);
    I1 = 2*pi*a(1)^2*a(3)/bb1^1.5*(acos(bb2)-bb2*(1-bb2^2)^0.5);
    I3 = 4*pi-2*I1;
    I13 = (I3-I1)/(3*bb1);
    I11 = pi/a(1)^2-3*I13/4;
    I12 = I11/3;
    I33 = 4*pi/(3*a(3)^2)-2*I13;
    oblate = [I11 I12 I13 I1;I12 I11 I13 I1;I13 I13 I33 I3];
end

function prolate = prolate(a)
    %Prolate
    bb1 = a(1)^2-a(3)^2;
    bb2 = a(1)/a(3);
    I2 = 2*pi*a(1)*a(3)^2/bb1^1.5*(bb2*(bb2^2-1)^0.5-acosh(bb2));
    I1 = 4*pi-2*I2;
    I12 = (I2-I1)/(3*bb1);
    I11 = 4*pi/(3*a(1)^2)-2*I12;
    I22 = pi/a(2)^2-3*I12/4;
    I23 = I22/3;
    prolate = [I11 I12 I12 I1;I12 I22 I23 I2;I12 I23 I22 I2];
end

function sphere = sphere(a)
    %Sphere
    I1 = 4*pi/3;
    I11 = 4*pi/(5*a(1)^2);
    I12 = I11/3;
    sphere = [I11 I12 I12 I1;I12 I11 I12 I1;I12 I12 I11 I1];
end

function Phi = Phi(a)
    %Phi
    theta = asin((1-(a(3)/a(1))^2)^0.5);
    k = ((a(1)^2-a(2)^2)/(a(1)^2-a(3)^2))^0.5;
    fun1 = @(w) 1./(1-k^2.*sin(w).^2).^0.5; 
    f = integral(fun1,0,theta);
    fun2 = @(w) (1-k^2.*sin(w).^2).^0.5;
    e = integral(fun2,0,theta);
    Phi = [f;e]; 
end
