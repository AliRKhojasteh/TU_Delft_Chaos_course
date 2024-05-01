function [La v] = iterate_coupledmap(ax,ay,epsi,Nit,x,y)


%%
lap = 0;
lu = -ax*log(ax) - (1-ax)*log(1-ax);

for i = 1:Nit-1
    
    if(x<=ax)
        fx = x/ax;
    else
        fx = (1-x)/(1-ax);
    end

    if(y<=ay)
        fy = y/ay;
    else
        fy = (1-y)/(1-ay);
    end
    
    
    x = (1-epsi)*fx + epsi*fy;
    y = epsi*fx + (1-epsi)*fy;

%% compute Lyapunov exponent statistics    
    v(i) = (x-y)/2;
    u = (x+y)/2;
    if(u<=ax)
        fd = 1/ax;
    else
        fd=1/(1-ax);
    end
    
    lap = lap + log(fd);
    
    La(i) = 1/i *lap - lu;
end


end