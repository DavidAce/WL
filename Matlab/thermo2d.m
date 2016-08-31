
function [ c,u,T, dos ] = thermo2d( dosEM,E,M,N )
%THERMO Summary of this function goes here
%   Detailed explanation goes here
T       = linspace(0,6,200)';
C		= zeros(length(T),1);
Z		= zeros(length(T),1);
u		= zeros(length(T),1);
dos     = zeros(length(E),1);
%%Subtract smallest value from dosE first


%Find first lambda to get g(E)
for j = 1:length(E)
    lambda = max(dosEM(j,:));
	dos(j) = lambda + log(nansum(exp(dosEM(j,:)-lambda))  );
end
for i = 1:length(T)
    t = T(i);
   %Find lambda
	lambda   = max(dos - E/t);
	DosExp   = exp(dos - E/t - lambda);
	Z(i)	 = nansum(DosExp);
	uAvg     = nansum(E.*DosExp)/Z(i);
	u(i)	 = uAvg/N;
	eSqAvg	 = nansum((E.^2).*DosExp)/Z(i);
	C(i) = (eSqAvg - uAvg*uAvg) /t^2 ;
end
c		= C/N;
end

