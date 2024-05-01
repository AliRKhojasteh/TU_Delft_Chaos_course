% template program using ODE45
clear all
close all

ga  = 0.2;
Gam = 1.2;
Om  = 0.5;
  
c   = 0.79;
tend = 100000;

% this is just the two coupled oscillators (1).  They are already made a little different
% (their gamma's differ by 1e-10).  To get the Lyapunov exponent of the single system, 
% you will have to add the single system, and its linearized version (2). To get the 
% Lyapunov exponent of the difference system, you will have to add the linearized difference
% system.  Just arrange it as in the assignment.
couplode = @(t,y) [ y(2); 
                    -(ga      )*y(2) - sin(y(1)) + Gam*cos(Om*t) + c*(sin(y(1)) - sin(y(3)))/2;
                    y(4);
                    -(ga+1e-10)*y(4) - sin(y(3)) + Gam*cos(Om*t) + c*(sin(y(3)) - sin(y(1)))/2]

                    % ode45 runs all the way to tend.  To get statistics of  finite-time 
                    % Lyapunov exponents you have to cut the interval [0 tend] in
                    % pieces.  For example 1e4 pieces of length 1e3.  After each piece, you
                    % set the initial condition to the last value.  From the normalization
                    % of the linearized vectors, you get the \Lambda_T.
                    
                    
                    % [t,y] = ode45(couplode, [0 tend], [0.1;0.1;0.1;0.1]); 
                    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
                    [t, y] = ode45(couplode, [0, tend], [0.1; 0.1; 0.1; 0.1], options);



                    % ode45 internally doeas as many time steps it needs, times are in t, solution in y.  'end' is the index of
                    % the final time, so y(end,:) is the final result.
                    
% just doing the differences for the histograms.
 
dif = sqrt( (y(:,1)-y(:,3)).^2 + (y(:,2)-y(:,4)).^2 );



% Plot the differences
figure;
plot(t, dif);
xlabel('Time');
ylabel('Differences');
