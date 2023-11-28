% program : quasip.m
%

% function B=quasip_map_0_picture_1()
    close all
    clear all

    global h alpha beta gamma A B Omega Picture Map
    
    animsteps = 32;
   
    
% Choose first between map and ode 
fprintf('Choose Map=1 for the map and Map=0 for the differential equation \n')
Map = 0

if(Map==0)   
    alpha     = 1;
    beta      = 1.576; 
    gamma     = 1; 
    A         = 1.4; 
    B         = 1; 
    Omega     = 1.76; 
    Picture   = 1; 
end

    fprintf('Press stop button to quit this run\n');

    h = 2*pi/(animsteps*Omega);         %step size
     
    xc = [0.3*2*pi 0.3 0];
    yy = xc(1)/(2*pi);
    tn = yy;
    rn =  0.3;
    t=0;
   % warning if Map is on in combination with picture 3
   if(Map==1 && Picture==3)
      fprintf('The option Map is not compatible with Picture option 3, Picture option is set to 1 \n');
      Picture = 1;
   end
%
H = init_figure();
    
nit = 20;
n = 1;

 while (n<nit)
        
    if(Map==0)
            b1 =  sin(xc(1));
            b2 = -cos(xc(1));
            B_coef(n,1)=b1;
            B_coef(n,2)=b2;

            % plot the current orientaion of the pendulum
            to = tn;
            xx = 0;
            h1=subplot(1,2,1);
            
           for i = 1:animsteps
                t = t + h;
                xc = Runge(xc);
                
               % plot pendulum
                cla(h1);
                b1 =  sin(xc(1));
                b2 = -cos(xc(1));
                plot(0, 0, '+', 'MarkerSize',10);
                plot([0 b1],[0 b2], 'b-');
                plot(b1, b2, 'r.', 'MarkerSize',25);
                pause(0.01)          
        end
        
             tn = xc(1)/(2*pi);  
             rn = xc(2);
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    if(Picture==1)
        % plot 1d map
        yy = mod(tn,1);
        xx = mod(to,1);
        if(Map==0) subplot(1,2,2);  end
        plot(xx,yy,'.k');
        pause(0.01)
    end  
   
    n = n+1;
 end
    

 error('you pressed stop :( ')

 % end



function H = init_figure()
    global Picture Map
    
    clf;
    H = uicontrol('Style', 'PushButton', ...
                    'String', 'Stop', ...
                    'Callback',@stopf);
                
if(Map==0)    
    
    subplot(1,2,1)
    axis([-1.5 1.5 -1.5 1.5]);
    axis('square');
    hold on
    xlabel('x');
    ylabel('y');
    title('Pendanim');

    subplot(1,2,2)
end
    axis([0 1 0 1])
    axis('square');
    hold on
    
    if(Picture==1)
    xlabel('$\theta_n$','interpreter','latex');
    ylabel('$\theta_{n+1}$','interpreter','latex');
    elseif(Picture==2)
    xlabel('$\theta_n$','interpreter','latex');
    ylabel('$\dot{\theta_n}$','interpreter','latex');
    if(Map==0)
    axis([0 1 -3 3])
    else
    axis([0 1 -0.5 0.5])
    end
    elseif(Picture==3)
    xlabel('t mod (2\pi / \Omega)');   
    ylabel('\theta mod 2\pi');    
    end
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.5]);
    commandwindow();  % Open Command Window, or select it if already open

end


function f = equations(x)
    global alpha beta gamma A B Omega
    f = zeros(3,1);
    f(1) = x(2);
    f(2) = 1/alpha * (-beta*x(2) - gamma*sin(x(1)) + A + B*cos(Omega*x(3)));
    f(3) = 1;
end


function xc = Runge(xc)
   global h
   
   x  = zeros(size(xc));
   c1 = zeros(size(xc));
   c2 = zeros(size(xc));
   c3 = zeros(size(xc));
   c4 = zeros(size(xc));

   n = length(xc);
   for i = 1:n; x(i) = xc(i); end
   f = equations(x);
   for i = 1:n; c1(i) = h*f(i); end

   for i = 1:n; x(i) = xc(i) + c1(i)/2; end
   f = equations(x);
   for i = 1:n; c2(i) = h*f(i); end

   for i = 1:n;  x(i) = xc(i) + c2(i)/2; end
   f = equations(x);
   for i = 1:n;  c3(i) = h*f(i); end

   for i = 1:n;  x(i) = xc(i) + c3(i); end
   f = equations(x);
   for i = 1:n;  c4(i) = h*f(i); end
   
   for i = 1:n
       xc(i) = xc(i) + (c1(i) + 2*c2(i) + 2*c3(i) + c4(i))/6;
   end
end   

function stopf(ObjectH, EventData)
delete(ObjectH);
end

