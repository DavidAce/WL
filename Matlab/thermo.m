function [ c,T ] = thermo( dos,E,N )
%THERMO Summary of this function goes here
%   Detailed explanation goes here
T = linspace(0,5,100);
C		= zeros(1,length(T));
Z		= zeros(1,length(T));
e		= zeros(1,length(T));

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

