function [ u ] = u_ising( T, L ) %
%Take a real number T, integer L and give U.
h = 1e-4;
fph = logZ_ising(T+h,L);
f   = logZ_ising(T,L);
u   = -2*L^2*sign(T) + T.^2.*(fph - f )./h;

% u = -2*L^2*sign(T) + T.^2.*gradient(logZ_ising(T,L),T);
u = u/L^2;
% u(isnan(u)) = 0;
end

