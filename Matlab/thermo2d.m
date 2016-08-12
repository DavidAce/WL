
function [ c,T ] = thermo2d( dosEM,E,M,N )
%THERMO Summary of this function goes here
%   Detailed explanation goes here
T = linspace(0,6,500);
C		= zeros(1,length(T));
Z		= zeros(1,length(T));
e		= zeros(1,length(T));
dos    = zeros(1,length(E));
%Find first lambda to get g(E)
for j = 1:length(E)
	lambda = 0;
	for l = 1:length(M)
		if dosEM(j,l) > lambda
			lambda = dosEM(j,l);
		end
	end
	dos(j) = lambda + log( sum(exp(dosEM(j,:)-lambda))  );
end
dos = dos';

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

