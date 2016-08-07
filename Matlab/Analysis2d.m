close all
clear all
for i = 0:3
    dos = importdata(['../outdata/dos' num2str(i) '.dat']);
    %dos(dos < 100) = NaN;

    E = importdata(['../outdata/E' num2str(i) '.dat']);
    M = importdata(['../outdata/M' num2str(i) '.dat']);
    N = 10^2;
    [c,T] = thermo2d(dos,E,M,N);
    figure(1)
    subplot(1,2,1)
    h = mesh(M,E, dos), hold on;
	dos(dos == 0) = NaN;
	dos(dos < 100) = NaN;
	set(h, 'zdata', dos);
	%axis vis3d
	ylabel('E');
	xlabel('M');
	zlabel('log(g(E,M))');
    axis([-100 100 -200 200])

    
%     surf(E,M,dos),hold on;

    subplot(1,2,2)
    plot(T,c),hold on;
end
    figure(5)
    h = mesh(M,E, dos);
	set(h, 'zdata', dos);
	axis vis3d
	ylabel('E');
	xlabel('M');
	zlabel('log(g(E,M))');