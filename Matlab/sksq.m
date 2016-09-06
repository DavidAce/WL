function [ val ] = sksq( alph ,beta ,k,m)
%CKSQ Summary of this function goes here
%   Detailed explanation goes here
sum = 0;
for j = 0:floor(m/2)
    sum = sum + factorial(m)/factorial(2*j)/factorial(m-2*j) * (alph(k)^2 - beta^2)^j * alph(k)^(m-2*j) ;
end

val = 1/2^(m-1)*(sum - beta^m);

end

