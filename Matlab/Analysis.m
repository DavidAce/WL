close all
clear all;
L = 6;
N = L^2;
cores = 4;
reps = 2;
for n = 1:reps

%      for i = 0:cores-1
%          
%         dos_sub = importdata(['../outdata/' num2str(n-1) '/dos' num2str(i) '.dat']);
%         %dos(dos < 5000) = 0;
%         E_sub = importdata(['../outdata/' num2str(n-1) '/E' num2str(i) '.dat']);
%         [c_sub,u_sub,T] = thermo(dos_sub,E_sub,N);
% 
%         figure(1);
%         set(gca,'ColorOrderIndex',i+1)
%         plot(E_sub,dos_sub, '-o'),hold all;
% 
%         figure(2);
%         set(gca,'ColorOrderIndex',i+1)
%         plot(T,c_sub),hold all;
% 
% %         figure(3);
% %         set(gca,'ColorOrderIndex',i+1)
% %         plot(T,u_sub),hold all;
%      end
     
    dos  = importdata(['../outdata/' num2str(n-1) '/dos.dat']);
    E    = importdata(['../outdata/' num2str(n-1) '/E.dat']);
    T    = importdata(['../outdata/' num2str(n-1) '/T.dat']);
    c    = importdata(['../outdata/' num2str(n-1) '/c.dat']);
    u    = importdata(['../outdata/' num2str(n-1) '/u.dat']);
    s    = importdata(['../outdata/' num2str(n-1) '/s.dat']);
    %f    = importdata(['../outdata/' num2str(n-1) '/f.dat']);

    %[c,u,T] = thermo(dos,E,N);
    cIsing = c_ising(T,L);
    uIsing = u_ising(T,L);
    cError = (cIsing - c)./max(c,cIsing);
    uError = (uIsing - u)./max(u,uIsing);


    figure(4)
    subplot(1,2,1)
    plot(T,c ,'DisplayName', ['c WL ' num2str(n)]),hold all;
    subplot(1,2,2)
    plot(T,cError ,'DisplayName', ['c error ' num2str(n)]),hold all;


    figure(5)
    subplot(1,2,1)
    plot(T,u ,'DisplayName', ['u WL ' num2str(n)]),hold all;
    subplot(1,2,2)
    plot(T,uError ,'DisplayName', ['u error' num2str(n)]),hold all;
    
    
    
    figure(6)
    plot(T,s ,'DisplayName', ['s WL ' num2str(n)]),hold all;
  
    figure(7)
%    plot(T,f ,'DisplayName', ['f WL ' num2str(n)]),hold all;
  
    

    [TIsing,eIsing,dosIsing] = dos_ising(E,L);
    dosIsing = dosIsing - (max(dosIsing) - max(dos));
    figure(8)
    subplot(1,2,1)
    plot(E/N,dos, '-o','DisplayName', ['dos WL ' num2str(n)]),hold all;

    subplot(1,2,2)
    dosError = (abs((dosIsing - dos)./(dosIsing)));
    semilogy(E/N,dosError),hold all;
end


figure(4)
subplot(1,2,1)
plot(T,cIsing, 'DisplayName','c Ising'),hold all;
xlabel('T')
ylabel('c(T)')
legend('show')

figure(5)
subplot(1,2,1)
plot(T,uIsing, 'DisplayName','u Ising'),hold all;
xlabel('T')
ylabel('u(T)')
legend('show')

figure(8)
subplot(1,2,1)
plot(eIsing, dosIsing, '-o', 'DisplayName','dos Ising');
legend('show')


