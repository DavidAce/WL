function [ c ] = c_ising( T, L ) %
%C_ISING Summary of this function goes here
%   Detailed explanation goes here
h = (1e-3);
fp2h = logZ_ising(T+2*h,L);
fph  = logZ_ising(T+h,L);
f    = logZ_ising(T,L);
c = 2*T.*(fph - f)/h/L^2 + T.^2.*(fp2h - 2*fph + f)/h^2/L^2;
% c = gradient(T.^2.*gradient(logZ_ising(T,L),T), T) / L^2 ;
% c = gradient(T.^2.*gradient(logZ_ising(T,L),T), T) / L^2 ;

end

