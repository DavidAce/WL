close all
clear all;
L = 10;
N = L^2;
dos1D = load(['../outdata/final/dos1D.dat']);
E     = load(['../outdata/final/E.dat']);
M     = load(['../outdata/final/M.dat']);
T     = load(['../outdata/final/T.dat']);
c     = load(['../outdata/final/c.dat']);
m     = load(['../outdata/final/m.dat']);
u     = load(['../outdata/final/u.dat']);
s     = load(['../outdata/final/s.dat']);
f     = load(['../outdata/final/f.dat']);
x     = load(['../outdata/final/x.dat']);
c_peak= load(['../outdata/final/c_peak.dat']);
x_peak= load(['../outdata/final/x_peak.dat']);
Tc_F  = load(['../outdata/final/Tc_F.dat']);
dos   = load(['../outdata/final/dos.dat']);
D     = load(['../outdata/final/D.dat']);
F     = load(['../outdata/final/F.dat']);
P     = load(['../outdata/final/P.dat']);


c_err     = load(['../outdata/final/c_err.dat']);
m_err     = load(['../outdata/final/m_err.dat']);
u_err     = load(['../outdata/final/u_err.dat']);
s_err     = load(['../outdata/final/s_err.dat']);
f_err     = load(['../outdata/final/f_err.dat']);
x_err     = load(['../outdata/final/x_err.dat']);
dos1D_err = load(['../outdata/final/dos1D_err.dat']);
c_peak_err= load(['../outdata/final/c_peak_err.dat']);
x_peak_err= load(['../outdata/final/x_peak_err.dat']);
Tc_F_err  = load(['../outdata/final/Tc_F_err.dat']);
dos_err   = load(['../outdata/final/dos_err.dat']);
D_err     = load(['../outdata/final/D_err.dat']);
F_err     = load(['../outdata/final/F_err.dat']);
P_err     = load(['../outdata/final/P_err.dat']);


%[c,u,T] = thermo(dos,E,N);
cIsing = c_ising(T,L);
cError = (cIsing - c)./max(c,cIsing);
uIsing = u_ising(T,L);
uError = (uIsing - u)./max(u,uIsing);

figure(1)
subplot(1,2,1)
shadedErrorBar(T,u,u_err,'-',1),hold all;
plot(T,uIsing);
title('Internal energy')
legend('u WL', 'u Ising')
xlabel('T');
ylabel('u(T)')

subplot(1,2,2)
shadedErrorBar(T,c, c_err, '-', 1),hold all;
errorbare('d', c_peak(1),c_peak(2), c_peak_err(1), c_peak_err(2))
plot(T,cIsing);
title('Specific heat')
legend('c WL', 'c Ising')
xlabel('T');
ylabel('c(T)')
% subplot(1,2,2)
% plot(T,cError ,'DisplayName', ['Ising Rel. Error ']),hold all;
% plot(T,c_err ,'DisplayName', ['Standard Error'])
% title('Error')
% xlabel('T');
% ylabel('Error')


figure(2)
subplot(1,2,1)
shadedErrorBar(T,s, s_err, '-',1),hold all;
title('Entropy')
legend('s WL')
xlabel('T');
ylabel('s(T)')
subplot(1,2,2)
shadedErrorBar(T,f ,f_err, '-',1),hold all;
title('Free Energy')
legend('f WL')
xlabel('T');
ylabel('f(T)')
%%
figure(3)
subplot(1,2,1)
shadedErrorBar(T,m ,m_err, '-',1),hold all;
title('Magnetization')
legend('m WL')
xlabel('T');
ylabel('m(T)')
subplot(1,2,2)
shadedErrorBar(T,x ,x_err, '-',1),hold all;
errorbare('d', x_peak(1),x_peak(2), x_peak_err(1), x_peak_err(2))
title('Susceptibility')
legend('x WL')
xlabel('T');
ylabel('x(T)')

figure(4)
[TIsing,eIsing,dosIsing] = dos_ising(E,L);
dosIsing = dosIsing - (max(dosIsing) - max(dos1D));
shadedErrorBar(E/N,dos1D,dos1D_err, '-o', 1),hold all;
plot(E/N,dosIsing);


figure(5)
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
figure(6)
subplot(1,2,1)
h = mesh(M/N,T, P);
hold on;
set(h, 'zdata', P);
axis vis3d
xlabel('M');
ylabel('T');
zlabel('log(g(E,M))');

subplot(1,2,2)
h = mesh(M/N,T, P_err);
hold on;
set(h, 'zdata', P_err);
axis vis3d
xlabel('M');
ylabel('T');
zlabel('log(g(E,M))');

%%
figure(7)
Tc_idx = nearestpoint( 2.6976205200244232074, T);
Tc_near= Tc_idx-20:10:Tc_idx+20;
shadedErrorBar(M/N, (F(Tc_idx,:)),F_err(Tc_idx,:)),hold all;
for i = 1:length(Tc_near)
plot(M/N,(F(Tc_near(i),:))),hold on
end
hold on;
xlabel('M');
ylabel('\Delta f(\beta,m)');
%%


figure(8)
Tc_idx = nearestpoint(c_peak(1), T);  %% Is the best canonical distribution really at this temperature?
shadedErrorBar(E, D(Tc_idx, :), D_err(Tc_idx,:), '-o', 1)

