clear all
close all
T = 2.26;
K = 1/T;
L = 10;
n = L;
m = L;
x = exp(-2*K);

beta =          2*x*(1-x^2)
alph = @(k)     (1+x^2)^2 - beta*cos(pi*k/n)
c0   =          (1-x)^m + x^m*(1+x)^m
s0   =          (1-x)^m - x^m*(1+x)^m
cn   =          (1+x)^m + x^m*(1-x)^m
sn   =          (1+x)^m - x^m*(1-x)^m
% 
% cksq_j = @(j) factorial(m)/factorial(2*j)/factorial(m-2*j) * (alph(k)^2 - beta^2)^j * alph(k)^(m-2*j); 
% cksq = zeros(floor(m/2) + 1);
% for j = 0:floor(m/2)
%     cksq(j) = cksq_j(j);
% end
% cksq = 1/2^(m-1)*(cksq + beta^m);
% 
% cksq_j = @(j) factorial(m)/factorial(2*j)/factorial(m-2*j) * (alph(k)^2 - beta^2)^j * alph(k)^(m-2*j); 
% cksq = zeros(floor(m/2) + 1);
% for j = 0:floor(m/2)
%     cksq(j) = cksq_j(j);
% end
% cksq = 1/2^(m-1)*(cksq + beta^m);



prod = 1;
for k = 0:n/2-1
    prod = prod * cksq(alph,beta,2*k+1,m);
end
Z1 = 0.5 * prod

prod = 1;
for k = 0:n/2-1
    prod = prod * sksq(alph,beta,2*k+1,m);
end
Z2 = 0.5 * prod

prod = 1;
for k = 1:n/2-1
    prod = prod * cksq(alph, beta,2*k,m);
end
Z3 = 0.5 * c0*cn * prod;

prod = 1;
for k = 1:n/2-1
    prod = prod * sksq(alph, beta,2*k,m);
end
Z4 = 0.5 * s0*sn * prod
syms k
cksq_exp = expand( cksq(alph, beta,2*k,m) );
sksq_exp = expand( sksq(alph, beta,2*k,m) );
log(cksq_exp) 


Z = exp(2*m*n*K) * ( Z1 + Z2 + Z3 + Z4)