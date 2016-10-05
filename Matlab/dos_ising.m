function [ T, e,S ] = dos_ising(E,L)
Tmax = 200;
Tmin = 1e-6;
% low = log(1e-6);
% high= log(1e4);
% d = (high - low)/n;
% T1 = exp(low:d:high)';


e = E(E <= 0)/L^2;
T = zeros(length(e),1);
parfor i = 1:length(e)
    T(i) = fminbnd(@(t) abs(u_ising(t,L) - e(i))  , 0 , Tmax);
end
e = [e; E(E>0)/L^2];
T = [T; -flip(T(1:end-1))];


f = -(2*L^2*sign(T) + T.*logZ_ising(T,L))/L^2;
S = (e - f)./T;
S = S*L^2;
% [T e f S logZ_ising(T,L)]
end

