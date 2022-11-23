function [vis] = eff_viscosity1_sI(A, Q, n, t, sI)
% find the effective viscosity from a known flow law (Ranalli,1987. P76-p78)
%   Based on flow law (1)
% input: flow law parameters
%         A:      preexponential parameter
%         Q:      activation energy (J/mol)
%         t:      temperature (Celcius)
%         n:      stress exponent
%         sI:     stress invariant
% output: vis:    effective viscosity

    T    = t +273.15; %K 
    R    = 8.314;  %J/mol/K
    A1   = 3^((n+1)/2)/2 * A;
    vis  = exp(Q/R/T)/2/A1*sI^(1-n); 
end