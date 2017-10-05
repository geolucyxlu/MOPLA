function q = Q(ang)
% Q.m
% 
% Transformation matrix Q from three spherical angles
%
% Input:  ang, three spherical angles defined in Jiang(2007a), in radian, 
%              3*1 matrix; 

% Output: q,   3*3 matrix
%--------------------------------------------------------------------------
    a = [sin(ang(2))*cos(ang(1)),sin(ang(2))*sin(ang(1)),cos(ang(2))];

    if(ang(2) == pi/2)
        b = [-sin(ang(1))*sin(ang(3)),cos(ang(1))*sin(ang(3)),cos(ang(3))];
    else
        b = [cos(atan(tan(ang(2))*cos(ang(1)-ang(3))))*cos(ang(3)),...
            cos(atan(tan(ang(2))*cos(ang(1)-ang(3))))*sin(ang(3)),...
             -sin(atan(tan(ang(2))*cos(ang(1)-ang(3))))];
    end
    c = cross(a,b);
    q = [a;b;c];
end