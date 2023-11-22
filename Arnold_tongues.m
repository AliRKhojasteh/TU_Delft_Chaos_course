%% Arnold tongues

clear all
close all


Niter = 1000;
t0 = 0;

Win = [0 1/2 3/5 8/13 5/8 2/3 1 Inf];
Wc(1:7)=0;
epsilon = 5e-04;
%% loop over omega and K starting from 0

K = [0:0.01:1.1];
lk = length(K);
Omega = [0:0.0001:1];
lo    = length(Omega);



for k = 1:lk
    n = 1;
    for o = 1:lo
        
        kf = K(end-k+1);
        om = Omega(o);

       t = t0;
       
       for s = 1:Niter
           t = t + om - kf/(2*pi) * sin(2*pi*t);
       end
       
       W = (t - t0)/s;
       
       if( abs(W-Win(n)) < epsilon )
           Kb(k,n) = kf;
           Ob(k,n) = om;
           Wp = Win(n);
           Wc(n) = 1;
           n = n+1;
       end
       
       if(Wc(n-1)==1)
           if(abs(W-Wp) > epsilon)
           Ke(k,n-1) = kf;
           Oe(k,n-1) = om;   
           Wc(n-1)=0;
           end
       end
        
    end
end


%%

figure
    plot(Ob,Kb,'.k');
    hold on
    plot(Oe,Ke,'.k');
