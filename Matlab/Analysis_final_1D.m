close all
clear all;
L = 16;
N = L^2;
dos1D= importdata(['../outdata/final/dos1D.dat']);
E    = importdata(['../outdata/final/E.dat']);
M    = importdata(['../outdata/final/M.dat']);
T    = importdata(['../outdata/final/T.dat']);
c    = importdata(['../outdata/final/c.dat']);
u    = importdata(['../outdata/final/u.dat']);
s    = importdata(['../outdata/final/s.dat']);
f    = importdata(['../outdata/final/f.dat']);
x    = importdata(['../outdata/final/x.dat']);
c_peak= importdata(['../outdata/final/c_peak.dat']);
dos   = importdata(['../outdata/final/dos.dat']);
D     = importdata(['../outdata/final/D.dat']);
F     = importdata(['../outdata/final/F.dat']);


c_err    = importdata(['../outdata/final/c_err.dat']);
u_err    = importdata(['../outdata/final/u_err.dat']);
s_err    = importdata(['../outdata/final/s_err.dat']);
f_err    = importdata(['../outdata/final/f_err.dat']);
x_err    = importdata(['../outdata/final/x_err.dat']);
dos1D_err = importdata(['../outdata/final/dos1D_err.dat']);
c_peak_err= importdata(['../outdata/final/c_peak_err.dat']);
dos_err   = importdata(['../outdata/final/dos_err.dat']);
D_err     = importdata(['../outdata/final/D_err.dat']);
F_err     = importdata(['../outdata/final/F_err.dat']);


%[c,u,T] = thermo(dos,E,N);
cIsing = c_ising(T,L);
cError = (cIsing - c)./max(c,cIsing);
uIsing = u_ising(T,L);
uError = (uIsing - u)./max(u,uIsing);

figure(1)
subplot(1,2,1)
shadedErrorBar(T,c, c_err, '-', 1),hold all;
errorbare('d', c_peak(1),c_peak(2), c_peak_err(1), c_peak_err(2))
plot(T,cIsing);
title('Specific heat')
legend('c WL', 'c Ising')
xlabel('T');
ylabel('c(T)')
subplot(1,2,2)
plot(T,cError ,'DisplayName', ['c error ']),hold all;
title('Relative Error')
xlabel('T');
ylabel('c error')

figure(2)
subplot(1,2,1)
shadedErrorBar(T,u,u_err,'-',1),hold all;
plot(T,uIsing);
title('Internal energy')
legend('u WL', 'u Ising')
xlabel('T');
ylabel('u(T)')
subplot(1,2,2)
plot(T,uError ,'DisplayName', ['u error']),hold all;
title('Relative Error')
xlabel('T');
ylabel('u error')

figure(3)
shadedErrorBar(T,s, s_err, '-',1),hold all;
title('Entropy')
legend('s WL')
xlabel('T');
ylabel('s(T)')


figure(4)
shadedErrorBar(T,f ,f_err, '-',1),hold all;
title('Free Energy')
legend('f WL')
xlabel('T');
ylabel('f(T)')

figure(5)
shadedErrorBar(T,x ,x_err, '-',1),hold all;
title('Susceptibility')
legend('x WL')
xlabel('T');
ylabel('x(T)')

% [TIsing,eIsing,dosIsing] = dos_ising(E,L);
% dosIsing = dosIsing - (max(dosIsing) - max(dos));
figure(6)
shadedErrorBar(E/N,dos1D,dos1D_err, '-o', 1),hold all;
% plot(E/N,dosIsing);
% subplot(1,2,2)
% dosError = (abs((dosIsing - dos)./(dosIsing)));
% semilogy(E/N,dosError),hold all;

%%
figure(7)
subplot(1,2,1)
h = mesh(M,E, dos);
hold on;
set(h, 'zdata', dos);
axis vis3d
ylabel('E');
xlabel('M');
zlabel('log(g(E,M))');
subplot(1,2,2)
h = mesh(M,E, dos_err);
hold on;
set(h, 'zdata', dos_err);
axis vis3d
ylabel('E');
xlabel('M');
zlabel('standard error');
%%

figure(8)
Tc_idx = nearestpoint(c_peak(1), T);
shadedErrorBar(E, D(Tc_idx, :), D_err(Tc_idx,:), '-o', 1)

