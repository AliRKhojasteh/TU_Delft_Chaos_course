%% Chaos week 6:    Synchronization

clear all
close all
%% Iterate coupled map

epscr = 0.2285;
ax = 0.7;
ay = 0.7+1e-9;

Nit = 1e6;
x = 0.5;
y = 0.5;

epsi = epscr-0.03;
[Llow vlow] = iterate_coupledmap(ax,ay,epsi,Nit,x,y);
epsi = epscr;
[Lcri vcri] = iterate_coupledmap(ax,ay,epsi,Nit,x,y);
epsi = epscr+0.03;
[Lup vup]  = iterate_coupledmap(ax,ay,epsi,Nit,x,y);
%%

% figure
% plot(x)
% hold on
% plot(y)

figure
plot(vlow)

z = log(abs(vlow/2));

[yn xn] = hist(z,50);
yn = yn/trapz(xn,yn);
figure
plot(xn,yn,'-ob')

[ma am] = max(yn);
xf = xn(1:am-1);
yf = yn(1:am-1);

%% fit histogram
xf(yf==0)=[];
yf(yf==0)=[];

figure
semilogy(xn,yn,'ob')

p = polyfit(xf(3:end),log(yf(3:end)),1)

hold on
semilogy(xf,exp(p(2))*exp(p(1)*xf),'-r')

%%

[yl xl] = hist(log(abs(vlow)),100);
yl = yl/trapz(xl,yl);
[yc xc] = hist(log(abs(vcri)),100);
yc = yc/trapz(xc,yc);
[yu xu] = hist(log(abs(vup)),100);
yu = yu/trapz(xu,yu);

%%
h=figure
semilogy(xl,yl,'-ob')
hold on
semilogy(xc,yc,'-ok')
hold on
semilogy(xu,yu,'-or')

ylabel P(z)
xlabel z

lg=legend('\epsilon_{cr}-0.03','\epsilon_{cr}','\epsilon_{cr}+0.03');
set(lg,'location','northwest')
legend boxoff

%print(h,'-dpng','Pz.png')