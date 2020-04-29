function [s] = SnP_poisson(a,nu)

% Eshelby Tensor Calculation for Interior points using elliptical integrals
% Note: p is not considered in this code.
    %
    s = zeros(3,3,3,3);
   % p = zeros(size(s));
    a = sort(a,'descend');
    j = J(a);
    %h = j(1:3,4);

    for lambda = 1:3
       for eta = 1:3
           if lambda==eta
               s(lambda,lambda,eta,eta) = (3/(8*pi*(1-nu)))*a(eta)^2*j(eta,eta) + ((1-2*nu)/(8*pi*(1-nu)))*j(eta,4);
           else
               s(lambda,lambda,eta,eta) = (1/(8*pi*(1-nu)))*a(eta)^2*j(lambda,eta) - ((1-2*nu)/(8*pi*(1-nu)))*j(lambda,4);
               s(lambda,eta,lambda,eta) = (1/(16*pi*(1-nu)))*(a(lambda)^2+a(eta)^2)*j(lambda,eta) + ((1-2*nu)/(16*pi*(1-nu)))*(j(lambda,4)+j(eta,4)) ;
               s(lambda,eta,eta,lambda) = s(lambda,eta,lambda,eta);
               
           end
           
       end
    end

    
% Note: All code based on Mura's and Qu and Cherkaoui's formalism
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
 %}
 %a=[5;3;1];
 %k= triaxial(a);
function triax = triaxial(a)
    % Triaxial
    
    phi= Phi(a);
    V = a(1)*a(2)*a(3);
    bb = [1/(a(1)^2-a(3)^2)^0.5; 1/(a(2)^2-a(3)^2)^0.5; 1/(a(1)^2-a(2)^2)^0.5];
    I1 = 4*pi*V*bb(1)*bb(3)^2*(phi(1)-phi(2));
    I3 = 4*pi*V*bb(1)*bb(2)^2*(a(2)*(a(1)^2-a(3)^2)^0.5/(a(1)*a(3))-phi(2));
    I2 = 4*pi-I1-I3;
    %% update: this part was modified to follow Mura's formulation
    % Date modified: 30/7/2019  by Ankit
    I12 = (I2-I1)*bb(3)^2; %/3;  
    I13 = (I3-I1)*bb(1)^2 ;%/3;
    I23 = (I3-I2)*bb(2)^2; %/3;
    I11 = (4*pi/(a(1)^2)-I12-I13)/3;  % here division by 3. Compare with SnP function based on Jiang's formulation for incompressible materials.
    I22 = (4*pi/(a(2)^2)-I12-I23)/3;
    I33 = (4*pi/(a(3)^2)-I13-I23)/3;
    triax = [I11 I12 I13 I1;I12 I22 I23 I2;I13 I23 I33 I3];
end
%
function oblate = oblate(a)
    % Oblate
    bb1 = a(1)^2-a(3)^2;
    bb2 = a(3)/a(1);
    I1 = 2*pi*a(1)^2*a(3)/bb1^1.5*(acos(bb2)-bb2*(1-bb2^2)^0.5);
    I3 = 4*pi-2*I1;
    %% update: this part was modified to follow Mura's formulation
    % Date modified: 30/7/2019  by Ankit
    I13 = (I3-I1)/(bb1);
    I11 = pi/a(1)^2-I13/4;
    I12 = I11;
    I33 = 4*pi/(3*a(3)^2)-2*I13/3;
    oblate = [I11 I12 I13 I1;I12 I11 I13 I1;I13 I13 I33 I3];
end

function prolate = prolate(a)
    %Prolate
    bb1 = a(1)^2-a(3)^2;
    bb2 = a(1)/a(3);
    I2 = 2*pi*a(1)*a(3)^2/bb1^1.5*(bb2*(bb2^2-1)^0.5-acosh(bb2));
    %% update: this part was modified to follow Mura's formulation
    % Date modified: 30/7/2019  by Ankit
    I1 = 4*pi-2*I2;
    I12 = (I2-I1)/(bb1);
    I11 = 4*pi/(3*a(1)^2)-2*I12/3;
    I22 = pi/a(2)^2-I12/4;
    I23 = I22;
    prolate = [I11 I12 I12 I1;I12 I22 I23 I2;I12 I23 I22 I2];
end

function sphere = sphere(a)
    %Sphere
    I1 = 4*pi/3;
    I11 = 4*pi/(5*a(1)^2);
    I12 = I11;
    sphere = [I11 I12 I12 I1;I12 I11 I12 I1;I12 I12 I11 I1];
end
%}
function Ph = Phi(a)
    %Phi
    theta = asin((1-(a(3)/a(1))^2)^0.5);
    k = ((a(1)^2-a(2)^2)/(a(1)^2-a(3)^2))^0.5;
    fun1 = @(w) 1./(1-k^2.*sin(w).^2).^0.5; 
    f = integral(fun1,0,theta);   % F(theta,k)in Jiang notes
    fun2 = @(w) (1-k^2.*sin(w).^2).^0.5;    
    e = integral(fun2,0,theta);   %E(theta,k) in Jiang notes
    Ph = [f;e]; 
end
end
