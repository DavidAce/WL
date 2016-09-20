clear all
close all

E     = importdata(['../outdata/final/E.dat']);
M     = importdata(['../outdata/final/M.dat']);
T     = importdata(['../outdata/final/T.dat']);
dos   = importdata(['../outdata/final/dos.dat']);


%Transform into 1d dos in M-space
lambda = max(max(dos));
loggM = log( lambda + nansum( exp(dos - lambda)))';
% lambdaM = max(dos)';
% for i = 1:length(M)
%     loggM(i) = log( lambdaM(i) + nansum(dos(:,i) - lambdaM(i)) );
% end

beta = 1./T;
betaE = T*E';

figure(1)
plot(M,loggM)


beta = 1./T;
betaE = T*E';

mexp =  nansum( M * exp( loggM - lambda - betaE)  )

