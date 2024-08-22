function lambda = solve_eq(abc, xp)

% Calculate a positive root of the equation
% x^2/(a^2+lambda)+ y^2/(b^2+lambda)+ z^2/(c^2+lambda)=1

% Input: abc are three semi-axes of the inclusion, 1*3 matrix;
%        xp are three coordinates of the points(n) we considered, defined in
%           C' system, 3*n matrix;

% Output: lambda is the maximum roots with corresponding xp, 1*n matrix

    abcs = abc.^2;   %square of a,b,c
    xyzs = xp.^2;   %square of x,y,z
    
    as = abcs(1);
    bs = abcs(2);
    cs = abcs(3);
    
    xs = xyzs(1,:);
    ys = xyzs(2,:);
    zs = xyzs(3,:);
    
    m  = as+bs+cs-(xs+ys+zs);
    n  = as*cs+bs*cs+as*bs-xs*(bs+cs)-ys*(as+cs)-zs*(as+bs);
    p  = as*bs*cs-(xs*bs*cs+ys*as*cs+zs*as*bs);
    tm = m.*n/6-m.^3/27-p/2;
    tn = n/3-m.^2/9;
    delt = tm.^2+tn.^3;
    
    r(1,:) = -m/3+(tm+delt.^0.5).^(1/3)+(tm-delt.^0.5).^(1/3);
    r(2,:) = -m/3+(-1+3^0.5j)*(tm+delt.^0.5).^(1/3)/2+(-1-3^0.5j)*(tm-delt.^0.5).^(1/3)/2;
    r(3,:) = -m/3+(-1-3^0.5j)*(tm+delt.^0.5).^(1/3)/2+(-1+3^0.5j)*(tm-delt.^0.5).^(1/3)/2;
    
    [~,n]  = size(xp);
    lambda = zeros(1,n);
    
    for i=1:n
             
        lambda(1,i) = real(max(real(r(:,i))));
      
    end
   
end