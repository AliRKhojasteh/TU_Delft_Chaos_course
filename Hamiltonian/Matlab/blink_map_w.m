
function blink_map()  % it is a function because we (could) callback
    
    clear all;
    xm = 15; ym = 15; % window size
      
    d_vortices = 4 * 3 * 72/25.4; % picture only, look in blink_w.m 

    % location of vortices
    xvort1 = -5; yvort1 = 0;
    xvort2 = 5;  yvort2 = 0;    

    T = 5 ; Omega = 10 ; 
    fprintf('\nDefault T, Omega %f %f\n', T, Omega);
    T = input('Give T = '); 
    if( isempty(T) ) T = 5; end
    Omega = input('give Omega = ');
    if( isempty(Omega) ) Omega = 10; end;
         
    Nsect = 5000;
    
    % arrays of stroboscopic points, Matlab has no way to tell it how
    % large, and that it should be double ... sigh
    % x_sect = linspace(0,1,Nsect) ; y_sect = linspace(0,1,Nsect);
    
    integration_steps = 2000; % so many integration steps in T
    h = T / integration_steps; % step size
  
    H = init_figure_blink(xm, ym); % make little GUI (with stop button, which is not needed)
    
    % plot the two vortex cores in blue   
    plot(xvort1, yvort1, 'b.', 'MarkerSize',d_vortices);
    plot(xvort2, yvort2, 'b.', 'MarkerSize',d_vortices);
   
    % You may want to make a loop around this, for example y_init = linspace(-10, 10, 8)
    % and thwn show the attractor for each y_init
    
    x_init = 0.0;  y_init = 2.0; % initial condition
    
    while get(H, 'Userdata') == 1 % stop button has not been pressed

    [x_init, y_init, button] = ginput(1); % mouse input on left click

    if(button==1) % initial point
        fprintf('Initial point: %6.3f, %6.3f\n', x_init, y_init);
    end
    
    isect = 1;  
  
    % this will also check for stop button (not needed, program is fast
    % enough)
    x = x_init; y = y_init;
    
  while((isect<=Nsect)&&(get(H,'Userdata') == 1))           
        
    if(mod(isect,2)==0) % if isect is even, the right (blue) vortex is on
        xc = xvort2 ; yc = yvort2;
    else % it is the left one
        xc = xvort1 ; yc = yvort1;
    end
    
    % integrate ! You could integrate two trajectories
    for istep = 1 : integration_steps-1 
        % simple forward Euler integration; there is an issue with accuracy
        % so you may want something more advanced
        x0 = x ; y0 = y ;
        r1sq = (x0 - xc)^2 + (y0 - yc)^2;
        x = x0 - Omega * (y0 - yc) / r1sq * h;
        y = y0 + Omega * (x0 - xc) / r1sq * h;
    end
          
    xsect(isect) = x ; ysect(isect) = y ;
    
    isect = isect + 1;

  end % section loop

% with these section coordinates you can do what you want
% now just plot them here

plot(xsect, ysect, '.r', 'Markersize', 3);

end

end


function H = init_figure_blink(xm, ym) 

    clf; % clear figure window
    % set up small user interface: just a stop button in the figure window 
    H = uicontrol('Style', 'PushButton', ...
                    'String', 'Stop', ...
                    'Callback',@stopf, 'Userdata', 1);
    
    axis([-xm xm -ym ym]);
    axis('square');
    title('Blinking point vortices');
    hold on;
    
    % size of window and where
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.7 0.7]);
    commandwindow();  % Open Command Window, or select it if already open
end

function stopf(ObjectH, EventData)
set(ObjectH, 'Userdata', 0);  % set Userdata 0 to stop
end

