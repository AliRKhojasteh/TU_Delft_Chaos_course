% must make this a function because of (elementary) GUI

function standard()
    clear all
      
    K = input('Give K = '); 
    if( isempty(K) ) K = 4; end

    niter = 5000;
    fprintf('Left mouse button gives %5d iterates, right is new window\n', niter);
    fprintf('Press stop button to quit altogether\n');

theta_min = 0 ; theta_max = 2*pi ; p_min = 0 ; p_max = 2*pi 

H = init_figure_standard(theta_min, theta_max, p_min, p_max);
    
twopi = 2*pi;

while get(H, 'Userdata') == 1 % stop button has not been pressed

[theta_0, p_0, button] = ginput(1); % mouse input on left click

if(button==1) % initial point
    fprintf('Initial point: %6.3f, %6.3f\n', theta_0, p_0);
end
if(button == 3) % new window
    theta_min = theta_0; 
    p_min = p_0;
    [theta_0, p_0, button] = ginput(1);
    theta_max = theta_0;
    p_max = p_0;
    
    % flip them in proper order, otherwise plot yells
    if(theta_min>theta_max) temp = theta_max ; theta_max = theta_min ; theta_min = temp; end
    if(p_min>p_max) temp = p_max ; p_max = p_min ; p_min = temp; end
       
    H = init_figure_standard(theta_min, theta_max, p_min, p_max);
    % replot what you have, but since this is zoomed, you will need 
    % much more
    plot(theta, p,'.r', 'Markersize', 3);
end          
n = 1;
p = zeros(1, n);
theta = zeros(1, n);

while (n <= niter)

    theta_1 = theta_0 + p_0;
    theta_1 = mod(theta_1, twopi);
    
    p_1 = p_0 + K * sin(theta_1);
    p_1 = mod(p_1, twopi);
    
    p(n) = p_1 ; theta(n) = theta_1; % store and plot later
    
    theta_0 = theta_1;
    p_0 = p_1;
    
    n = n+1;  
end
    plot(theta, p,'.r', 'Markersize', 3);
    % you now have an array of niter iterates, do with it what you want
end
  
 end


 function H = init_figure_standard(theta_min,theta_max,p_min,p_max)
       
    clf; % clear figure window
    % set up small user interface: just a stop button in the figure window 
    H = uicontrol('Style', 'PushButton', ...
                    'String', 'Stop', ...
                    'Callback',@stopf, 'Userdata', 1);
                
    axis([theta_min theta_max p_min p_max])
    axis('square');
    hold on
    
    xlabel('$\theta_n$','interpreter','latex');
    ylabel('$p_{n+1}$','interpreter','latex');
       
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.7 0.7]);
    commandwindow();  % Open Command Window, or select it if already open

end


function stopf(ObjectH, EventData)
set(ObjectH, 'Userdata', 0);
end

