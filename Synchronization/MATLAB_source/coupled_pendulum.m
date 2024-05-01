%% differential equations

clear all
close all

%% set parameters
global ga Gam Om c h
    
ga  = 0.2;
Gam = 1.2;
Om  = 0.5;

c   = 0.79;
tm = [0.1 0.2 0];
ts = [0.1 0.2 0];
Nit = 100000;
h = 0.001;

%%

for k = 1:Nit
    
    tm = Runge_master(tm);  
    ts = Runge_slave(ts,tm);
    
    et(k) = log(sqrt( (ts(1)-tm(1))^2 + (ts(2)-tm(2))^2 ));
    v(k) = sqrt( (ts(1)-tm(1))^2 + (ts(2)-tm(2))^2 );
    tmm(1:3,k) = tm;
    tss(1:3,k) = ts;
    
end

%%

[yn xn] = hist(et,25);

figure
semilogy(xn,yn,'-ob')

p = polyfit(xn(4:12),log(yn(4:12)),1)

hold on
semilogy(xn,exp(p(2))*exp(xn*p(1)),'-r')


