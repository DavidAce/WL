function [ c,e,T ] = thermo( dos,E,N )
%THERMO Summary of this function goes here
%   Detailed explanation goes here
T = linspace(1e-6,6,500)';
C		= zeros(length(T),1);
Z		= zeros(length(T),1);
e		= zeros(length(T),1);

for i = 1:length(T)
    t = T(i);
   %Find lambda
	lambda = max(dos - E/t);
	DosExp = exp(dos - E/t - lambda);
	Z(i)	 = sum(DosExp);
	eAvg     = sum(E.*DosExp)/Z(i);
	e(i)	 = eAvg/N;
	eSqAvg	 = sum((E.^2).*DosExp)/Z(i);
	C(i) = (eSqAvg - eAvg*eAvg) /t^2 ;
end
c		= C/N;

end

