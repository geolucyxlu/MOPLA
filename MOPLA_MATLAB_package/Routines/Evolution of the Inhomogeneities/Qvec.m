function q = Qvec(ang)
% Qvec.m
% Vectorized version of Q function 
%
%--------------------------------------------------------------------------
    [~,n] = size(ang);
    a(1,:) = sin(ang(2,:)).*cos(ang(1,:));
    a(2,:) = sin(ang(2,:)).*sin(ang(1,:));
    a(3,:) = cos(ang(2,:));

    id    = find(ang(2,:) == pi/2); 
    othid = find(ang(2,:) ~= pi/2);
    
    b       = zeros(3,n);
    b(1,id) = -sin(ang(1,id)).*sin(ang(3,id));
    b(2,id) = cos(ang(1,id)).*sin(ang(3,id));
    b(3,id) = cos(ang(3,id));
    b(1,othid) = cos(atan(tan(ang(2,othid)).*cos(ang(1,othid)-ang(3,othid)))).*cos(ang(3,othid));
    b(2,othid) = cos(atan(tan(ang(2,othid)).*cos(ang(1,othid)-ang(3,othid)))).*sin(ang(3,othid));
    b(3,othid) = -sin(atan(tan(ang(2,othid)).*cos(ang(1,othid)-ang(3,othid)))); 

    c = cross(a,b);
    q = zeros(3,3,n);

    q(1,:,:) = a;
    q(2,:,:) = b;
    q(3,:,:) = c;

end