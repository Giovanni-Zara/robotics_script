function [T0N, rotated_jacobians] = direct_kinematics(N, DHTABLE, J, verbose)
    % Input:
    % N: Number of joints
    % DHTABLE: DH parameters table
    % J: Jacobian matrix (MxN) expressed in the base frame (M can be any number of rows)

    % Definire variabili simboliche con nomi che non confliggono con nomi predefiniti
    syms alpha_sym d_sym a_sym theta_sym;

    % Costruire la matrice di trasformazione Denavit-Hartenberg generica
    TDH = [ cos(theta_sym) -sin(theta_sym)*cos(alpha_sym)  sin(theta_sym)*sin(alpha_sym) a_sym*cos(theta_sym);
            sin(theta_sym)  cos(theta_sym)*cos(alpha_sym) -cos(theta_sym)*sin(alpha_sym) a_sym*sin(theta_sym);
              0             sin(alpha_sym)                cos(alpha_sym)                d_sym;
              0               0                              0                          1];

    % Creare una cella per le matrici di trasformazione
    A = cell(1, N);

    % Sostituire i valori nella matrice DH generica per ogni giunto
    for i = 1:N
        alpha_sym = DHTABLE(i, 1);
        a_sym = DHTABLE(i, 2);
        d_sym = DHTABLE(i, 3);
        theta_sym = DHTABLE(i, 4);
        A{i} = subs(TDH, {'alpha_sym', 'a_sym', 'd_sym', 'theta_sym'}, {alpha_sym, a_sym, d_sym, theta_sym});
    end

    T = eye(4);

    % Initialize cell array to store rotation matrices and rotated Jacobians
    rotation_matrices = cell(1, N);
    rotated_jacobians = cell(1, N);

    for i = 1:N 
        T = T * A{i};
        T = simplify(T);

        % Extract the rotation matrix R_{0i} from T_{0i}
        R = T(1:3, 1:3);
        rotation_matrices{i} = R;
        if verbose == true
            disp(['0_R_', num2str(i), ' ALREADY TRANSPOSED:']);
            disp(rotation_matrices{i}')
        end

        % Determine the number of rows in the Jacobian
        M = size(J, 1);

        % Rotate the Jacobian to the i-th frame
        if M >= 3
            % Rotate the linear part (first 3 rows)
            J_linear = J(1:3, :);
            J_linear_rotated = R' * J_linear;
            J_rotated = J_linear_rotated;

            if M > 3
                % Rotate the angular part (rows 4 to 6)
                J_angular = J(4:6, :);
                J_angular_rotated = R' * J_angular;
                J_rotated = [J_linear_rotated; J_angular_rotated];
            end
        else
            % If the Jacobian has fewer than 3 rows, assume it's only the linear part
            J_rotated = R' * J;
        end

        rotated_jacobians{i} = J_rotated;
        
        disp(['Rotated Jacobian in frame ', num2str(i), ':']);
        if isa(J_rotated, 'sym')
            disp(simplify(J_rotated));
        else
            disp(J_rotated);
        end
        disp("----------------------------------------------------------------------")

    end

    % Output the final transformation matrix T0N
    T0N = T;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

N = 3;
syms q1 q2 q3 a1 d1 a2 a3 real;
DHTABLE =  [ pi/2 a1 d1 q1;
                0 a2  0 q2;
             pi/2 a3  0 q3
            ];

% Construct the Jacobian matrix J_L(q)
J = [-sin(q1)*(a1 + a3*cos(q2 + q3) + a2*cos(q2)), -cos(q1)*(a3*sin(q2 + q3) + a2*sin(q2)), -a3*sin(q2 + q3)*cos(q1);
 cos(q1)*(a1 + a3*cos(q2 + q3) + a2*cos(q2)), -sin(q1)*(a3*sin(q2 + q3) + a2*sin(q2)), -a3*sin(q2 + q3)*sin(q1);
                                          0,            a3*cos(q2 + q3) + a2*cos(q2),          a3*cos(q2 + q3)];

verbose = true;

[T0N, rotated_jacobians] = direct_kinematics(N, DHTABLE, J, verbose);
