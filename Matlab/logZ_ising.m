function [ logZ ] = logZ_ising( T_vec, L )
logZ = (zeros(length(T_vec),1));
for i = 1:length(T_vec)
T = abs(T_vec(i)); 
K = (1/T);
n = (L);
m = (L);
x = (exp(-2*K));

beta =          (2*x*(1-x^2));
alph = @(k)     ((1+x^2)^2 - beta*cos(pi*k/n));
c0   =          ((1-x)^m + x^m*(1+x)^m);
s0   =          ((1-x)^m - x^m*(1+x)^m);
cn   =          ((1+x)^m + x^m*(1-x)^m);
sn   =          ((1+x)^m - x^m*(1-x)^m);

prod = 1;
for k = 0:n/2-1
    prod = prod * cksq(alph,beta,2*k+1,m);
end
Z1 = 0.5 * prod;

prod = 1;
for k = 0:n/2-1
    prod = prod * sksq(alph,beta,2*k+1,m);
end
Z2 = 0.5 * prod;

prod = 1;
for k = 1:n/2-1
    prod = prod * cksq(alph, beta,2*k,m);
end
Z3 = 0.5 * c0*cn * prod;

prod = 1;
for k = 1:n/2-1
    prod = prod * sksq(alph, beta,2*k,m);
end
Z4 = 0.5 * s0*sn * prod;
logZ(i) =  log(Z1 + Z2 + Z3 + Z4);
end

end

