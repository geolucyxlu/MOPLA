function [S,PI] = SnPI_el(a,nu)
% Calculation of interior Eshelby tensors for isotropic elastic materials
% using elliptic functions.
% input: a: three semi axes of the ellipsoidal inclusion
%        nu: Poisson's ratio
% output: Eshleby tensors S and PI
% LL 10 Oct,2023

S = zeros(3,3,3,3);
PI = zeros(size(S));
JJ = J(a);
h = JJ(1:3,4);
q = 3/(8*pi*(1-nu));
r = (1-2*nu)/(8*pi*(1-nu));
for i = 1:3
   for j = 1:3
       if i ~= j
           S(i,j,i,j) = q/2*(a(i)^2+a(j)^2)*JJ(i,j) + r/2*(h(i)+h(j));
           S(i,j,j,i) = S(i,j,i,j);
           S(i,i,j,j) = q*a(j)^2*JJ(i,j) - r*h(i);

           PI(i,j,i,j) = (h(j)-h(i))/(8*pi);
           PI(i,j,j,i) = PI(i,j,i,j);
           PI(j,i,i,j) = -PI(i,j,i,j);
           PI(j,i,j,i) = -PI(i,j,i,j);
       else
           S(i,i,i,i)  = q*a(j)^2*JJ(i,j) + r*h(i);
           PI(i,i,i,i) = 0;
       end
   end
end

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
