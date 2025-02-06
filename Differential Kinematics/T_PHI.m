%% Main function
function T = compute_T(rotation_type, rotation_sequence)
    if strcmpi(rotation_type, 'RPY')
        T = compute_T_RPY_CASE(rotation_sequence);
    else
        T = compute_T_EULER_CASE(rotation_sequence);
    end
end

%% Function for Euler Case
function T_return = compute_T_EULER_CASE(rotation_sequence)
    disp("Rotation Type: EULER");
    disp(['Sequence: ', rotation_sequence]);
    R_cumulative = eye(3);  % Placeholder for the cumulative rotation matrix
    T_return = sym(zeros(3,3));

    % Iterate over the rotation sequence
    for i = 1:length(rotation_sequence)
        axis = rotation_sequence(i); % Current rotation axis
        syms init alp bet gam;
        sign = init;
        switch i
            case 1
                sign = alp;
            case 2
                sign = bet;
            case 3
                sign = gam;
        end

        switch axis
            case 'X'
                unit_vector = [1; 0; 0];
            case 'Y'
                unit_vector = [0; 1; 0];
            case 'Z'
                unit_vector = [0; 0; 1];
            otherwise
                error('Invalid axis in rotation sequence.');
        end
        
        column_vector = R_cumulative * unit_vector;
        T_return(:, i) = column_vector;
       
        %%rimuovere commento se si vuole vedere evoluzione matrice
        %disp(T_return);

        switch axis
            case 'X'
                R_current = [1, 0, 0; 0, cos(sign), -sin(sign); 0, sin(sign), cos(sign)];
            case 'Y'
                R_current = [cos(sign), 0, sin(sign); 0, 1, 0; -sin(sign), 0, cos(sign)];
            case 'Z'
                R_current = [cos(sign), -sin(sign), 0; sin(sign), cos(sign), 0; 0, 0, 1];
        end

        % Multiply the cumulative rotation matrix with the current rotation matrix
        R_cumulative = R_cumulative * R_current; 
    end
end

%% RPY CASE
function T_return = compute_T_RPY_CASE(rotation_sequence)
    disp("Rotation Type: RPY");
    rotation_sequence = rotation_sequence(end:-1:1);
    disp(['Reversed sequence: ', rotation_sequence]);
    % Placeholder for the cumulative rotation matrix
    R_cumulative = eye(3); 
    T_return = sym(zeros(3,3));

    % Iterate over the rotation sequence
    for i = 1:length(rotation_sequence)
        axis = rotation_sequence(i); % Current rotation axis
        syms init alp bet gam;
        sign = init;
        switch i
            case 1
                sign = gam;
            case 2
                sign = bet;
            case 3
                sign = alp;
            otherwise
                error('Invalid angle in rotation sequence.');
        end

        switch axis
            case 'X'
                unit_vector = [1; 0; 0];
            case 'Y'
                unit_vector = [0; 1; 0];
            case 'Z'
                unit_vector = [0; 0; 1];
            otherwise
                error('Invalid axis in rotation sequence.');
        end
        
        column_vector = R_cumulative * unit_vector; 
        T_return(:, 4-i) = column_vector;
                
        %%rimuovere commento se si vuole vedere evoluzione matrice
        %disp(T_return);

        switch axis
            case 'X'
                R_current = [1, 0, 0; 0, cos(sign), -sin(sign); 0, sin(sign), cos(sign)];
            case 'Y'
                R_current = [cos(sign), 0, sin(sign); 0, 1, 0; -sin(sign), 0, cos(sign)];
            case 'Z'
                R_current = [cos(sign), -sin(sign), 0; sin(sign), cos(sign), 0; 0, 0, 1];
        end
        
        % Multiply the cumulative rotation matrix with the current rotation matrix
        R_cumulative = R_cumulative * R_current; 
    end
end

%% Function to substitute the sign with the numerical value 
function matrix_substituted = substitute_params(matrix, alpha_val, beta_val, gamma_val)
    % Substitutes the given values in the symbolic matrix for alpha, beta, gamma
    syms alp bet gam real;
    if ~isempty(alpha_val)
        matrix = subs(matrix, alp, alpha_val);
    end
    if ~isempty(beta_val)
        matrix = subs(matrix, bet, beta_val);
    end
    if ~isempty(gamma_val)
        matrix = subs(matrix, gam, gamma_val);
    end
    % Check if the matrix is fully numeric
    if ~has(matrix, [alp, bet, gam])
        % Convert to double if no symbolic variables remain
        matrix_substituted = double(matrix);
    else
        % Keep the matrix as symbolic if symbolic variables remain
        disp('---Matrix NOT FULLY NUMERIC, still contains symbolic terms---');
        matrix_substituted = matrix;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% COMPUTING THE T MATRIX
type_rotation = 'RPY'; % 'RPY' or 'Euler'
% WARNING THE FOLLOWING SEQUENCE MUST BE EQUAL TO THE ONE IN THE TEXT
AXIS_sequence = 'XZY';  % use CapsLock!!
T_matrix = compute_T(type_rotation, AXIS_sequence);
fprintf("THE FOLLOWING MATRIX IS T(PHI), will be used to compute w = T(PHI) * PHI' \n\n")
disp(T_matrix);

% IN CASE IS NEEDED TO ASSIGN CERTAIN VALUES IN THE T MATRIX
alpha_val = [];  % Leave as [] if not initialized
beta_val  = []; % Leave as [] if not initialized
gamma_val = [];  % Leave as [] if not initialized

disp('Substituted matrix:');
matrix_substituted = substitute_params(T_matrix, alpha_val, beta_val, gamma_val);
disp(matrix_substituted);

