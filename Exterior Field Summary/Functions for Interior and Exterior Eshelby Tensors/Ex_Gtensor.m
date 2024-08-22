function G= Ex_Gtensor(a,x)
% Ex_Gtensor.m
% Exterior-point Eshelby tensor (G_bar) calculation based on Gao & Ma,2010
%
% input: a: 3 semi-axes of the Eshelby ellipsoid (3*1 vector).
%        x: coordinates of exterior points (3*N matrix).
% output: G
%--------------------------------------------------------------------------

lambda = solve_eq(a, x); %lambda: 1*N

q      = x(1,:).^2./(lambda+a(1)^2).^2 + x(2,:).^2./(lambda+a(2)^2).^2 +...
         x(3,:).^2./(lambda+a(3)^2).^2; %q: 1*N
q1     = x(1,:).^2./(lambda+a(1)^2).^3 + x(2,:).^2./(lambda+a(2)^2).^3 +...
         x(3,:).^2./(lambda+a(3)^2).^3; %q1: 1*N
omega  =[(lambda+a(1).^2).^-1; (lambda+a(2).^2).^-1; (lambda+a(3).^2).^-1];       
         %omega: 3*N

s1     = S1(lambda, a);                      %s1: 3*3*N
s2     = S2(lambda, a);                      %s2: 3*3*N
s3     = S3(lambda, a, q, omega);            %s3: 3*N
s4     = S4(lambda, a, q, omega);            %s4: 3*N 
s5     = S5(lambda, a, q, omega);            %s5: 3*N 
s7     = S7(lambda, a, q, q1, omega);        %s7: 3*3*3*3*N

[~,n]  = size(lambda);
G      = zeros(3,3,3,3,n);
delt   = eye(3);
for i=1:3
    for j=1:3
        for l=1:3
            for m=1:3
                t1 = s1(i,l,:)*delt(i,j)*delt(l,m) + s2(i,j,:)*...
                    (delt(i,l)*delt(j,m)+delt(i,m)*delt(j,l)); %1*1*N
                tt1= (squeeze(t1))';
                t2 = s3(i,:).*delt(i,j).*omega(l,:).*omega(m,:).*x(l,:)...
                    .*x(m,:);%1*N
                t3 = s4(l,:).*delt(l,m).*omega(i,:).*omega(j,:).*x(i,:)...
                    .*x(j,:);%1*N 
                t4 = s5(i,:).*(omega(j,:).*omega(m,:).*delt(i,l).*x(j,:)...
                    .*x(m,:) + omega(j,:).*omega(l,:).*delt(i,m).*x(j,:)...
                    .*x(l,:));%1*N
                t5 = s5(j,:).*(omega(i,:).*omega(m,:).*delt(j,l).*x(i,:)...
                    .*x(m,:) + omega(i,:).*omega(l,:).*delt(j,m).*x(i,:)...
                    .*x(l,:));%1*N
                t6 = (squeeze(s7(i,j,l,m,:)))'.*x(i,:).*x(j,:).*x(l,:)...
                    .*x(m,:);%1*N
                G(i,j,l,m,:) = tt1 + t2 + t3 + t4 + t5 + t6;
                G(j,i,l,m,:) =  G(i,j,l,m,:);
                G(i,j,m,l,:) =  G(i,j,l,m,:);
                G(j,i,m,l,:) =  G(i,j,l,m,:); 
            end
        end
    end
end
end

function J = J1(lambda, a)
%J: 3*N
[~,n]  = size(lambda);
u1     = [lambda+a(1).^2; lambda+a(2).^2; lambda+a(3).^2];
u2     = [a(1)^2-a(2)^2; a(1)^2-a(3)^2; a(2)^2-a(3)^2];
u3     = repmat(u2,1,n);
u      = cat(1,u1,u3);
theta  = asin((u(5,:)./u(1,:)).^0.5);  
k      = (u(4,:)./u(5,:)).^0.5;
f      = lellipf(theta, k, 1.0e-10);
e      = lellipe(theta, k, 1.0e-10);
aa     = (f-e).*2./u(4,:)./u(5,:).^0.5;
b      = 2.* u(5,:).^0.5./u(4,:)./u(6,:).*e - 2./u(4,:)./u(5,:).^0.5.*f - ...
         2.* u(3,:).^0.5./u(6,:)./(u(1,:).*u(2,:)).^0.5;
c      = 2./u(6,:).*(u(2,:)./u(1,:)./u(3,:)).^0.5 - 2./u(6,:)./u(5,:).^0.5.*e;
t      = cat(1,aa,b,c);
J      = 2*pi*a(1)*a(2)*a(3).*t;

end

function JJ = J2(lambda, a)
%JJ: 3*3*N
[~,n]  = size(lambda);
u1     = [lambda+a(1).^2; lambda+a(2).^2; lambda+a(3).^2];
u2     = [a(1)^2-a(2)^2; a(1)^2-a(3)^2; a(2)^2-a(3)^2];
u3     = repmat(u2,1,n);
u      = cat(1,u1,u3);
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
k      = J1(lambda, a)./(2*pi*a(1)*a(2)*a(3));
v      = zeros(3,3,n);
for i=1:3
    for j=1:3
        if i~=j
           v(i,j,:)= (k(i,:)-k(j,:))/(a(j)^2-a(i)^2);
        end
    end
end
v(1,1,:) = (-2/3)./u(1,:)./delta + ((2*a(1)^2-a(2)^2-a(3)^2)/3)...
           ./u(5,:)./u(4,:).*k(1,:) - (1/3)./u(4,:).*k(2,:) - ...
           (1/3)./u(5,:).*k(3,:);
v(2,2,:) = (-2/3)./u(2,:)./delta - ((2*a(2)^2-a(1)^2-a(3)^2)/3)...
           ./u(6,:)./u(4,:).*k(2,:) + (1/3)./u(4,:).*k(1,:) - ...
           (1/3)./u(6,:).*k(3,:);    
v(3,3,:) = (-2/3)./u(3,:)./delta + ((2*a(3)^2-a(1)^2-a(2)^2)/3)...
           ./u(5,:)./u(6,:).*k(3,:) + (1/3)./u(5,:).*k(1,:) + ...
           (1/3)./u(6,:).*k(2,:);
JJ      = (2*pi*a(1)*a(2)*a(3)).*v; 
end

function s1 = S1(lambda, a)
%s1: 3*3*N
[~,n]  = size(lambda);
J      = J1(lambda, a);
JJ     = J2(lambda, a);
s1     = zeros(3,3,n);
for i=1:3
    for j=1:3
       s1(i,j,:) = (1/4/pi).*(J(i,:)-J(j,:)+a(i)^2.*(squeeze(JJ(i,j,:)))');
    end
end
end

function s2 = S2(lambda, a)
%s2: 3*3*N
[~,n]  = size(lambda);
J      = J1(lambda, a);
JJ     = J2(lambda, a);
s2     = zeros(3,3,n);
for i=1:3
    for j=1:3
      s2(i,j,:) = (1/4/pi).*((J(i,:)-J(j,:))./2+a(i)^2.*(squeeze(JJ(i,j,:)))');
    end
end
end

function s3 = S3(lambda, a, q, omega)
%s3: 3*N
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
t      = (a(1)*a(2)*a(3)).*lambda./delta./q; %t:1*N
tt     = repmat(t,3,1); 
s3     = tt.*omega;
end

function s4 = S4(lambda, a, q, omega)
%s4: 3*N 
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
t      = (a(1)*a(2)*a(3))./delta./q;
t1     = repmat(t,3,1);
t2     = (repmat(lambda,3,1).*omega-1); 
s4     = t1.*t2;
end

function s5 = S5(lambda, a, q, omega)
%s5: 3*N 
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
t      = (a(1)*a(2)*a(3))./delta./q;
t1     = repmat(t,3,1);
t2     = (repmat(lambda,3,1).*omega-0.5);
s5     = t1.*t2;
end

function s7 = S7(lambda, a, q, q1, omega)
%s7: 3*3*3*3*N
[~,n]  = size(lambda);
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
theta  = sum(omega);
v      = zeros(3,3,n);
s7     = zeros(3,3,3,3,n);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                v(k,l,:)=(a(1)*a(2)*a(3))./delta./q.^2.*...
                          omega(i,:).*omega(j,:).*omega(k,:).*omega(l,:).*...
                          (2 - 2.*lambda.*(omega(i,:)+omega(j,:)+omega(k,:)...
                          +omega(l,:))- lambda.*theta + 4.*lambda.*q1./q);
            end
        end
        s7(i,j,:,:,:) = v;
    end
end
end
