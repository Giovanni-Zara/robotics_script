function [posLaw, velLaw, accelLaw] = derive_or_integrate(baseLaw, type)
    % This function computes the position, velocity, or acceleration laws
    % starting from one of them (position, velocity, or acceleration).
    %
    % Inputs:
    %   baseLaw - Anonymous function representing the base law (pos, vel, accel).
    %   type    - A string specifying the input type ('pos', 'vel', or 'accel').
    %
    % Outputs:
    %   posLaw  - Position law (if calculable, symbolic and anonymous).
    %   velLaw  - Velocity law (if calculable, symbolic and anonymous).
    %   accelLaw- Acceleration law (if calculable, symbolic and anonymous).

    syms t C1 C2 % Define time variable and constants of integration
    lawSym = baseLaw(t); % Convert anonymous function to symbolic

    switch type
        case 'pos'
            % Starting from position -> derive for velocity and acceleration
            posLaw = simplify(lawSym); % Position is given
            velLaw = simplify(diff(lawSym, t)); % Velocity = d(position)/dt
            accelLaw = simplify(diff(lawSym, t, 2)); % Acceleration = d^2(position)/dt^2
        case 'vel'
            % Starting from velocity -> integrate for position, derive for acceleration
            velLaw = simplify(lawSym); % Velocity is given
            accelLaw = simplify(diff(lawSym, t)); % Acceleration = d(velocity)/dt
            posLaw = simplify(int(lawSym, t) + C1); % Position = integral(velocity) + C1
        % old version
        % case 'accel'
        %     % Starting from acceleration -> integrate for velocity and position
        %     accelLaw = simplify(lawSym); % Acceleration is given
        %     velLaw = simplify(int(lawSym, t) + C1); % Velocity = integral(acceleration) + C1
        %     posLaw = simplify(int(int(lawSym, t) + C1, t) + C2); % Position = double integral(acceleration) + C1 + C2
        case 'accel'
            % Starting from acceleration -> integrate for velocity and position
            accelLaw = sym(baseLaw(t)); % Convert anonymous function to symbolic
            velLaw = simplify(int(accelLaw, t) + C1); % Velocity = integral(acceleration) + C1
            posLaw = simplify(int(velLaw, t) + C2); % Position = integral(velocity) + C2
        otherwise
            error('Invalid input type. Use "pos", "vel", or "accel".');
    end
    
    %% SE CI SONO DELLE VARIABILI NELL'EQUAZIONI (SU CUI NON DEVI DERIVARE O INTEGRARE DEVI INSERIRLE 
    %ESPLICITAMENTE(DOPO AVERLE DICHIARATE NELLA FUNZIONE) NELLE ULTIME TRE RIGHE DELLA FUNZIONE

    % Return both symbolic and anonymous functions
    if exist('posLaw', 'var'), posLaw = struct('symbolic', posLaw, 'function', matlabFunction(posLaw, 'Vars', {t, C1, C2})); end
    if exist('velLaw', 'var'), velLaw = struct('symbolic', velLaw, 'function', matlabFunction(velLaw, 'Vars', {t, C1})); end
    if exist('accelLaw', 'var'), accelLaw = struct('symbolic', accelLaw, 'function', matlabFunction(accelLaw, 'Vars', t)); end
end

%% SE CI SONO DELLE VARIABILI NELL'EQUAZIONI (SU CUI NON DEVI DERIVARE O INTEGRARE DEVI INSERIRLE 
%ESPLICITAMENTE(DOPO AVERLE DICHIARATE NELLA FUNZIONE) NELLE ULTIME TRE RIGHE DELLA FUNZIONE

%% STARTING FROM POSITION
% posLaw = @(t) 0.5 * 200 * t.^2;
% 
% % Compute the velocity and acceleration laws
% [~, velLaw, accelLaw] = derive_or_integrate(posLaw, 'pos');
% 
% % Test the output
% disp("---------------------------------------ATTENZIONE---------------------------------------");
% disp("RICODARSI EVENTUALI CONDIZIONI INIZIALI SE NON SI è AL PRIMO PERIODOOO")
% disp("----------------------------------------------------------------------------------------");
% disp('Velocity Law:');
% disp(velLaw);
% disp('Acceleration Law:');
% disp(accelLaw);
% 
% disp('Readable Laws:');
% fprintf('Velocity Law: v(t) = %s\n', char(velLaw.symbolic));
% fprintf('Acceleration Law: a(t) = %s\n', char(accelLaw.symbolic));

%% STARTING FROM VELOCITY
% Define the velocity law
% velLaw = @(t) 200 * t;
% 
% % Compute the position and acceleration laws
% [posLaw, ~, accelLaw] = derive_or_integrate(velLaw, 'vel');
% 
% % Display symbolic laws
% disp("---------------------------------------ATTENZIONE---------------------------------------");
% disp("RICODARSI EVENTUALI CONDIZIONI INIZIALI SE NON SI è AL PRIMO PERIODOOO")
% disp("----------------------------------------------------------------------------------------");
% disp('Symbolic Laws:');
% disp(['Position Law: ', char(posLaw.symbolic)]);
% disp(['Acceleration Law: ', char(accelLaw.symbolic)]);
% 
% % Evaluate the position law at t = 1, assuming C1 = 0
% t_val = 1;
% C1 = 0;
% pos_at_t1 = posLaw.function(t_val, C1);
% disp(['Position at t = 1 (C1 = 0): ', num2str(pos_at_t1)]);

%% STARTING FROM ACCELERATION
% Define the acceleration law
accelLaw = @(t) 200 * ones(size(t));

% Compute the position and velocity laws
[posLaw, velLaw, ~] = derive_or_integrate(accelLaw, 'accel');

disp("---------------------------------------ATTENZIONE---------------------------------------");
disp("RICODARSI EVENTUALI CONDIZIONI INIZIALI SE NON SI è AL PRIMO PERIODOOO")
disp("----------------------------------------------------------------------------------------");
% Test the output
disp('Position Law:');
disp(posLaw);

disp('Velocity Law:');
disp(velLaw);
