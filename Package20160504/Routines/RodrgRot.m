function z = RodrgRot(A)
% RodrgRot.m
% Rodrigues' Rotation Approximation
% Input:  A, 3*3 matrix;

% Output: z, 3*3 matrix; 
%--------------------------------------------------------------------------
    omega = norm(A)/2^0.5;
    if omega == 0
        z = eye(3);
    else
       omega1 = A/omega;
       z      = eye(3) + omega1 * sin(omega) + (1 - cos(omega)) * omega1^2;
    end
end