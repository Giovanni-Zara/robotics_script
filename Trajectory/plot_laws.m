function plot_general_joint_trajectory(laws, periods, toPlot, total_time)
    % Inputs:
    % laws      - A structure containing anonymous functions for acceleration, velocity, and position:
    %             laws(k).accel = @(t) ... (acceleration law for period k)
    %             laws(k).vel = @(t) ... (velocity law for period k)
    %             laws(k).pos = @(t) ... (position law for period k)
    % periods   - An array specifying the time intervals for each law: [t_start, t_end; ...]
    % toPlot    - Cell array specifying which plots to create, e.g., {'pos', 'vel', 'accel'}
    % total_time - The total expected time for the trajectory
    
    % Check input validity
    if length(laws) ~= size(periods, 1)
        disp(length(laws));
        disp(size(periods, 1));
        disp('------ERROREEEE-----')
        error('The number of laws must match the number of time intervals.');
    end
    
    % Check if the sum of period lengths matches total_time
    period_lengths = periods(:, 2) - periods(:, 1);
    if abs(sum(period_lengths) - total_time) > 1e-6
        fprintf('ERROREEEEE sum of period lengths (%.2f) does not match the total_time (%.2f).', sum(period_lengths), total_time);
        error('SOMMA PERIODI NON COINCIDE CON TOTAL TIME');
    end

    % Prepare combined time and trajectory data
    t_combined = [];
    accel_combined = [];
    vel_combined = [];
    pos_combined = [];

    % Iterate over each period and evaluate laws
    for k = 1:length(laws)
        % Extract time interval for this period
        t_start = periods(k, 1);
        t_end = periods(k, 2);
        t = linspace(t_start, t_end, 100); % Fine-grained time vector for plotting

        % Evaluate the laws
        accel = laws(k).accel(t); % Acceleration
        vel = laws(k).vel(t);     % Velocity
        pos = laws(k).pos(t);     % Position

        % Append results
        t_combined = [t_combined, t];
        accel_combined = [accel_combined, accel];
        vel_combined = [vel_combined, vel];
        pos_combined = [pos_combined, pos];
    end

    % Plot results based on user selection
    if any(strcmp(toPlot, 'pos'))
        plot_position(t_combined, pos_combined);
    end
    if any(strcmp(toPlot, 'vel'))
        plot_velocity(t_combined, vel_combined);
    end
    if any(strcmp(toPlot, 'accel'))
        plot_acceleration(t_combined, accel_combined);
    end
end

% Function to plot position
function plot_position(t, pos)
    figure;
    plot(t, pos, 'r', 'LineWidth', 1.5);
    xlabel('time [s]');
    ylabel('Position');
    title('Joint Trajectory (Position)');
    grid on;
end

% Function to plot velocity
function plot_velocity(t, vel)
    figure;
    plot(t, vel, 'b', 'LineWidth', 1.5);
    xlabel('time [s]');
    ylabel('Velocity');
    title('Joint Velocities');
    grid on;
end

% Function to plot acceleration
function plot_acceleration(t, accel)
    figure;
    plot(t, accel, 'g', 'LineWidth', 1.5);
    xlabel('time [s]');
    ylabel('Acceleration');
    title('Joint Accelerations');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TIPOLOGIE DI LAWS
% laws(1).accel = @(t) 200 * ones(size(t));                       % Constant acceleration
% laws(1).vel = @(t) 200 s* t;                                     % Linear velocity
% laws(1).pos = @(t) 0.5 * 200 * t.^2;                            % Quadratic position
% 
% laws(2).accel = @(t) -200 * ones(size(t));                      % Constant deceleration
% laws(2).vel = @(t) 200 * (1 - t);                               % Linear velocity (reverse)
% laws(2).pos = @(t) 100 - 0.5 * 200 * (1 - t).^2;                % Quadratic position (reverse)




%IMPORTANTE CANCELLARE SEMPREEEEE
clear laws periods; 

% Define laws for two time periods
laws(1).accel = @(t) 200 * ones(size(t));                       % Constant acceleration
laws(1).vel = @(t) 200 * t;                                     % Linear velocity
laws(1).pos = @(t) 0.5 * 200 * t.^2;                            % Quadratic position

laws(2).accel = @(t) -200 * ones(size(t));                      % Constant deceleration
laws(2).vel = @(t) 200 * (1 - t);                               % Linear velocity (reverse)
laws(2).pos = @(t) 100 - 0.5 * 200 * (1 - t).^2;                % Quadratic position (reverse)

% Define time intervals for the four laws
total_time = 2;
first_period = [0, 1];
second_period = [1, 2];

% Combine them into a matrix
periods = [first_period; second_period];

% Specify which plots to generate
toPlot = {'pos', 'vel', 'accel'}; % Choose any combination: 'pos', 'vel', 'accel'

% Plot the trajectory
plot_general_joint_trajectory(laws, periods, toPlot, total_time);
