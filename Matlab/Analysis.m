close all
for i = 0:3
    dos = importdata(['../outdata/dos' num2str(i) '.dat']);
    %dos(dos < 5000) = 0;
    E = importdata(['../outdata/E' num2str(i) '.dat']);
    N = 10^2;
    [c,T] = thermo(dos,E,N);

    figure(1)
    plot(E,(dos), '-o'),hold on;

    figure(2)
    plot(T,c),hold on;
end
