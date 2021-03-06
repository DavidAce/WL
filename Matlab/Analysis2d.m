close all
clear all
L = 12;
N = L^2;
for j = 0:4
    for i = 0:1:8
        dos = load(['../outdata/' num2str(j) '/dos' num2str(i) '.dat']);
        E   = load(['../outdata/' num2str(j) '/E'   num2str(i) '.dat']);
        M = load(  ['../outdata/' num2str(j) '/M' num2str(i) '.dat']);
        [c,u,T, dosE] = thermo2d(dos,E,M,N);
        figure(10+j);
        subplot(1,2,1);
        dos(dos==0) = nan;
        dos = dos - min(min(dos));
        h = mesh(M,E, dos); 
        hold on;

        set(h, 'zdata', dos);
        ylabel('E');
        xlabel('M');
        zlabel('log(g(E,M))');
        %axis([-100 100 -200 200])


    %     surf(E,M,dos),hold on;

        subplot(1,2,2)
        plot(T,c),hold on;
    end
end
%%
    figure(2);
    dos = importdata(['../outdata/0/dos.dat']);
    E = importdata(['../outdata/0/E.dat']);
    M = importdata(['../outdata/0/M.dat']);
    [c,u,T,dosE] = thermo2d(dos,E,M,N);
    subplot(1,2,1);
    h = mesh(M,E, dos);
    hold on;
 	set(h, 'zdata', dos);
	axis vis3d
	ylabel('E');
	xlabel('M');
	zlabel('log(g(E,M))');
    
    cIsing = c_ising(T,L);
    uIsing = u_ising(T,L);
    cError = ((cIsing - c));
    UError = ((uIsing - u));
    
    subplot(1,2,2);
    plot(T,c,'-o', 'DisplayName','c WL'),hold on;
    plot(T,cIsing,'-o', 'DisplayName','c EX');
    ylabel('c');
	xlabel('T');
    legend('show')
    figure(3);
    plot(T,u),hold on;
    plot(T,uIsing);
    ylabel('U');
	xlabel('T');
    figure(4);
    subplot(1,2,1);
    plot(T,UError);
    ylabel('U error');
	xlabel('T');
    subplot(1,2,2);
    plot(T,cError);
    ylabel('c error');
	xlabel('T');
    %%
    figure(5)
    [TIsing,eIsing,dosIsing] = dos_ising(E,L);
    %dosIsing = dosIsing + (max(dosE) - max(dosIsing));
    plot(N*eIsing,dosIsing),hold on 
    plot(E,dosE);
    
    figure(6)
    plot(E,interp1(N*eIsing,dosIsing,E) - dosE)
    