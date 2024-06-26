%% differential equations

clear all
close all

%% set parameters
global ga Gam Om c h
    
ga  = 0.2;
Gam = 1.2;
Om  = 0.5;

c   = 0.79;
tend = 100000;

cr = c;
%%

couplode = @(t,y) [ y(2); 
                    -ga*y(2) - sin(y(1)) + Gam*cos(Om*t) + c*(sin(y(1)) - sin(y(3)))/2;
                    y(4);
                    -(ga+1e-10)*y(4) - sin(y(3)) + (Gam)*cos(Om*t) + c*(sin(y(3)) - sin(y(1)))/2]

[t,y] = ode45(couplode, [0 tend], [0.1;0.1;0.1;0.1]); 
%%
%                 
% diffode = @(tt,d) [d(2);
%                    -ga*d(2) + 2*(1-cr)*d(1)*cos(y())]
%                
% 
% [tt,d] = ode45(diffode, [0 tend], [1e-8;0.1]);             
 %%
 
 dif = sqrt( (y(:,1)-y(:,3)).^2 + (y(:,2)-y(:,4)).^2 );
 eta = log(dif);
 
 %%
 figure
 plot(t,dif);
 ylabel('\eta(t)')
 xlabel t
 
 [yn xn] = hist(eta,50);
 yn = yn./trapz(xn,yn);
 
 figure
 semilogy(xn,yn,'ob')
 xlabel('\eta');
 ylabel('P(\eta)')
 
 %%
 figure
 plot(t,y(:,[1 3]))
 figure
 plot(t,y(:,[2 4])) 
 
 figure
 plot(tt,d)