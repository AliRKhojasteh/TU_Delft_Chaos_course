field_size = 0.5;
n_points = 100;

% Blue fluid
x1 = linspace(-field_size, field_size, n_points);
y1 = linspace(-field_size, field_size, n_points);



[X1, Y1] = meshgrid(x1, y1);
% [X2, Y2] = meshgrid(x2, y2);

t_final = 10000;
n_steps = 20000;
dt = t_final / n_steps;

% Vortex positions
v1_pos = [-1, 0];
v2_pos = [1, 0];

delta = .5;

% Frequency for switching vortices
vortex_frequency = 50;
orientation = 1;

% Video setup
video = VideoWriter('fluid_simulation.avi');
open(video);

for t = 1:n_steps
    % Switch vortices
    if mod(t, vortex_frequency) == 0
        orientation = orientation * -1;
    end
    
    % Update positions using velocity field
    if orientation == 1
        [u1, v1] = velocity_field(X1, Y1, v1_pos(1), v1_pos(2), orientation, delta);
        X1 = X1 + dt * u1;
        Y1 = Y1 + dt * v1;
        
        % [u2, v2] = velocity_field(X2, Y2, v1_pos(1), v1_pos(2), orientation, delta);
        % X2 = X2 + dt * u2;
        % Y2 = Y2 + dt * v2;
    elseif orientation == -1
        [u1, v1] = velocity_field(X1, Y1, v2_pos(1), v2_pos(2), orientation, delta);
        X1 = X1 + dt * u1;
        Y1 = Y1 + dt * v1;
        
        % [u2, v2] = velocity_field(X2, Y2, v2_pos(1), v2_pos(2), orientation, delta);
        % X2 = X2 + dt * u2;
        % Y2 = Y2 + dt * v2;
    end
    
    scatter(X1, Y1, 'red');
    % hold on;
    % scatter(X2, Y2, 'blue');
    hold off;

    xlim([-2, 2]);
    ylim([-2, 2]);
    % axis off;

    drawnow;

    % Capture the plot as an image
    frame = getframe(gcf);
    writeVideo(video, frame);

    clf; % Clear figure after capturing the frame
end

close(video);

% Function computing velocity field for the shear flow
function [u, v] = velocity_field(x, y, v_posx, v_posy, orientation, delta)
    u = (-orientation / (2 * pi)) * ((y - v_posy) ./ ((x - v_posx) .^ 2 + (y - v_posy) .^ 2 + delta ^ 2));
    v = (orientation / (2 * pi)) * ((x - v_posx) ./ ((x - v_posx) .^ 2 + (y - v_posy) .^ 2 + delta ^ 2));
end
