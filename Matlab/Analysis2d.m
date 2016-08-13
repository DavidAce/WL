close all
clear all
N = 10^2;

for i = 0:5
    dos = importdata(['../outdata/dos' num2str(i) '.dat']);
    dos(dos > 0) = dos(dos>0) + randi(100);

    E = importdata(['../outdata/E' num2str(i) '.dat']);
    M = importdata(['../outdata/M' num2str(i) '.dat']);
    [c,T] = thermo2d(dos,E,M,N);
    figure(1)
    subplot(1,2,1)
    h = mesh(M,E, dos), hold on;
	dos(dos == 0) = NaN;
	set(h, 'zdata', dos);
	%axis vis3d
	ylabel('E');
	xlabel('M');
	zlabel('log(g(E,M))');
    %axis([-100 100 -200 200])

    
%     surf(E,M,dos),hold on;

    subplot(1,2,2)
    plot(T,c),hold on;
end
    figure(2)
    dos = importdata(['../outdata/dos.dat']);
    E = importdata(['../outdata/E.dat']);
    M = importdata(['../outdata/M.dat']);
    [c,T] = thermo2d(dos,E,M,N);
    subplot(1,2,1)
    h = mesh(M,E, dos), hold on;
	dos(dos == 0) = NaN;
	set(h, 'zdata', dos);
	%axis vis3d
	ylabel('E');
	xlabel('M');
	zlabel('log(g(E,M))');
    subplot(1,2,2)
    plot(T,c),hold on;