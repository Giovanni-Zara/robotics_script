
%% Calcola la matrice di trasformazione DH, supportando parametri simbolici
% alpha: angolo di rotazione attorno all'asse x
% d: traslazione lungo l'asse z
% a: traslazione lungo l'asse x
% theta: angolo di rotazione attorno all'asse z
function A = dh_matrix(alpha, d, a, theta)
    A = [cos(theta), -sin(theta) * cos(alpha),  sin(theta) * sin(alpha), a * cos(theta);
         sin(theta),  cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta);
         0,           sin(alpha),              cos(alpha),              d;
         0,           0,                       0,                       1];
end

%% Funzione per moltiplicare una serie di matrici in ordine.
% varargin: serie di matrici da moltiplicare
%---------------------------------------------------------------------
% esempio utilizzo funzione
% syms a11 a12 a21 a22 b11 b12 b21 b22 c11 c12 c21 c22
% A = [a11, a12; a21, a22];
% B = [b11, b12; b21, b22];
% C = [c11, c12; c21, c22];
% result = multiply_matrices(A, B, C); ----> vanno inserite cosi
%---------------------------------------------------------------------
function result = multiply_matrices(varargin)
    if nargin < 1
        error('Devi fornire almeno una matrice come input.');
    end
    result = varargin{1};
    for i = 2:nargin
        result = result * varargin{i};
    end
    result = simplify(result);
end

%% Funzione che restituisce la matrice di rotazione ATTORNO A UNO DEI TRE ASSI e all'angolo specificato
% usi le tre matrici di rotazoni classiche per ruotare intorno a ASSI x, y, z
% axis: asse di rotazione COME STRINGA ('x', 'y', 'z')
% angle: angolo di rotazione SIA COME NUMERIC CHE SIMBOLICO
function R = rotation_matrix(axis, angle)
    % Check if the axis is valid
    if ~ismember(axis, {'x', 'y', 'z'})
        error('Asse non valido. Usa ''x'', ''y'', o ''z''.');
    end
  
    % Ensure the angle is symbolic if it is not numeric
    if ~isa(angle, 'sym')
        angle = sym(angle);
    end
    
    % Define the rotation matrix based on the specified axis
    switch axis
        case 'x'
            R = [1, 0, 0;
                 0, cos(angle), -sin(angle);
                 0, sin(angle), cos(angle)];
        case 'y'
            R = [cos(angle), 0, sin(angle);
                 0, 1, 0;
                 -sin(angle), 0, cos(angle)];
        case 'z'
            R = [cos(angle), -sin(angle), 0;
                 sin(angle), cos(angle), 0;
                 0, 0, 1];
        otherwise
            error('Asse non valido. Usa ''x'', ''y'', o ''z''.');
    end
end

%% Funzione che restituisce p_hom date le matrici di trasformazione DH
% varargin: serie di matrici di trasformazione DH (3x3) da moltiplicare
function result = compute_p_hom(varargin)
    col = [0; 0; 0; 1];
    result = multiply_matrices(varargin{:});
    result = result * col;
    result = simplify(result);
end

%% Funzione che estrae il vettore invariante rispetto alla matrice di rotazione
% R: matrice di rotazione
function result = not_rotated_vector(R)
    [V, D] = eig(R); % autovettori in V, autovalori in D
    index = diag(D)==1;
    result = V(:, index);
    result = vpa(result);
end

%% Funzione che normalizza un vettore
% vec: vettore da normalizzare
function result = normalize_vector(vec)
    result = vec/norm(vec);
    result = vpa(result);
end

%% Funzione che calcola il sin di theta dalla matrice di rotazione R_theta_r
% per il problema INVERSO 
% R: matrice di rotazione
function result = compute_sin(R)
    result = (1/2) * sqrt((R(1, 2) - R(2, 1))^2 + (R(1, 3) - R(3, 1))^2 + (R(2, 3) - R(3, 2))^2);
    if result == 0
        warning('WARNING: Sine value is zero. SINGULAR CASE.');
    end
end

%% Funzione che calcola il cos di theta dalla matrice di rotazione R_theta_r
% per il problema INVERSO 
% R: matrice di rotazione
function result = compute_cos(R)
    result = (1/2) * (R(1, 1) + R(2, 2) + R(3, 3) -1);
end

%% Funzione che calcola il vettore r data la matrice di rotazione R_theta_r e sin(theta)
% per il problema INVERSO 
% R: matrice di rotazione
% s: sin(theta)
function result = compute_r(R, s)
    arr = [R(3, 2) - R(2, 3);
           R(1, 3) - R(3, 1);
           R(2, 1) - R(1, 2)];
    result = (1/(2*s)) * arr;
end

%% Funzione che calcola la matrice di rotazione di angolo theta attorno al VETTORE v
% ---> è esattamente il direct problem scritto in una altra maniera
% intorno a VETTORE v (NON asse x, y, z)
% v: vettore attorno al quale ruotare
% theta: angolo di rotazione
function R = rotation_matrix_for_r(v, theta)
    if ~isa(theta, 'sym')
        theta = sym(theta);
    end

    v = v / norm(v);

    vx = v(1);
    vy = v(2);
    vz = v(3);
    
    R = [cos(theta) + vx^2 * (1 - cos(theta)),   vx * vy * (1 - cos(theta)) - vz * sin(theta),   vx * vz * (1 - cos(theta)) + vy * sin(theta);
         vy * vx * (1 - cos(theta)) + vz * sin(theta),   cos(theta) + vy^2 * (1 - cos(theta)),   vy * vz * (1 - cos(theta)) - vx * sin(theta);
         vz * vx * (1 - cos(theta)) - vy * sin(theta),   vz * vy * (1 - cos(theta)) + vx * sin(theta),   cos(theta) + vz^2 * (1 - cos(theta))];
    
    R = simplify(R);
end

%% Funzione che calcola la matrice di rotazione utilizzando gli angoli di roll, pitch e yaw
% a_1, a_2, a_3: assi di rotazione ('x', 'y', 'z')
% alpha_1, alpha_2, alpha_3: angoli di rotazione
function result = roll_pitch_yaw(a_1, a_2, a_3, alpha_1, alpha_2, alpha_3)
    if ~ismember(a_1, {'x', 'y', 'z'})
        error('Asse 1 non valido. Usa ''x'', ''y'', o ''z''.');
    end
    if ~ismember(a_2, {'x', 'y', 'z'})
        error('Asse 2 non valido. Usa ''x'', ''y'', o ''z''.');
    end
    if ~ismember(a_3, {'x', 'y', 'z'})
        error('Asse 3 non valido. Usa ''x'', ''y'', o ''z''.');
    end

    R_1 = rotation_matrix(a_1, alpha_1);
    R_2 = rotation_matrix(a_2, alpha_2);
    R_3 = rotation_matrix(a_3, alpha_3);

    result = multiply_matrices(R_3, R_2, R_1);
end

%% Funzione che calcola la matrice di rotazione utilizzando gli angoli di Eulero
% a_1, a_2, a_3: assi di rotazione ('x', 'y', 'z')
% alpha_1, alpha_2, alpha_3: angoli di rotazione
function result = euler_angles(a_1, a_2, a_3, alpha_1, alpha_2, alpha_3)
    if ~ismember(a_1, {'x', 'y', 'z'})
        error('Asse 1 non valido. Usa ''x'', ''y'', o ''z''.');
    end
    if ~ismember(a_2, {'x', 'y', 'z'})
        error('Asse 2 non valido. Usa ''x'', ''y'', o ''z''.');
    end
    if ~ismember(a_3, {'x', 'y', 'z'})
        error('Asse 3 non valido. Usa ''x'', ''y'', o ''z''.');
    end

    R_1 = rotation_matrix(a_1, alpha_1);
    R_2 = rotation_matrix(a_2, alpha_2);
    R_3 = rotation_matrix(a_3, alpha_3);

    result = multiply_matrices(R_1, R_2, R_3);
end

%% Funzione che calcola l'angolo theta 
% sin_val: valore di sin(theta)
% cos_val: valore di cos(theta)
function angle = custom_atan2(sin_val, cos_val)
    % Print the signs of the input parameters
    if sin_val > 0
        disp('The sine value is positive.');
    elseif sin_val < 0
        disp('The sine value is negative.');
    else
        disp('The sine value is zero.');
    end

    if cos_val > 0
        disp('The cosine value is positive.');
    elseif cos_val < 0
        disp('The cosine value is negative.');
    else
        disp('The cosine value is zero.');
    end

    % Calculate the angle using atan2
    angle = atan2(sin_val, cos_val);
end

%% Direct kinematics 
% N: number of joints
% DHTABLE: table containing the DH parameters 

% -------------------- esempio di utilizzo di questa funzione --------------------
% syms q_1 q_2 l1 l2
% N = 2
% DHTABLE =  [pi -l1 0 q_1; 
%            -pi/2 l2 0 q_2 ];
% [T0N, p, n, s, a] = direct_kinematics(N, DHTABLE)
%--------------------------------------------------------------------------------

function [T0N, p, n, s, a] = direct_kinematics(N, DHTABLE)
    
    disp(['Number of joints N=', num2str(N)])
    disp('DH table')

    % Definire variabili simboliche con nomi che non confliggono con nomi predefiniti
    syms alpha_sym d_sym a_sym theta_sym

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

    for i = 1:N 
        T = T * A{i};
        T = simplify(T);
    end

    % Output della matrice T0N
    disp('output O_T_N matrix')
    T0N = T;
    % Output della posizione
    disp('output ON position')
    p = T(1:3, 4);
    % Output degli assi
    disp('normal axis output (x)')
    n = T(1:3, 1);
    disp('stride axis output (y)')
    s = T(1:3, 2);
    disp('approach axis output (z)')
    a = T(1:3, 3);
end


%% Computes the rotation matrix given a vector r and an angle theta.
% r (vector): a 3-dimensional vector.
% theta: the rotation angle in radians, can be numeric or symbolic.
function R = solve_direct_problem(r, theta)

    % Ensure r is a column vector and normalize it if needed
    r = r(:);  % Convert to a column vector if not already

    % Calculate the norm of r and normalize if necessary
    norm_r = sqrt(r' * r);
    if norm_r ~= 1
        disp('Normalizing vector r');
        r = r / norm_r;  % Normalize the vector
    end

    % Define symbolic variables for theta if not numeric
    if ~isa(theta, 'sym')
        theta = sym(theta);
    end

    % Identity matrix
    identity_matrix = eye(3);
    % r * r^T (outer product)
    r_transposed = r * r';
    % Extract elements for the cross-product matrix
    r_x = r(1);
    r_y = r(2);
    r_z = r(3);

    % Skew-symmetric cross-product matrix
    s_of_r_matrix = [0, -r_z, r_y;
                     r_z, 0, -r_x;
                     -r_y, r_x, 0];
    % Compute the rotation matrix using Rodrigues' rotation formula
    R = r_transposed + (identity_matrix - r_transposed) * cos(theta) + s_of_r_matrix * sin(theta);
    R = simplify(R);
end

%% Check if the input matrix is a rotation matrix 
function is_rotation_matrix(R)
    % A rotation matrix should be orthonormal, meaning R' * R = I and det(R) = 1
    
    % Check if the matrix is square
    [rows, cols] = size(R);
    if rows ~= cols
        disp('La matrice non è quadrata.');
        disp('Risultato: Non è una matrice di rotazione.');
        return;
    end
    
    % Check if the matrix is orthonormal (R' * R should be close to I)
    I = eye(rows);
    R_transpose_R = R' * R;
    if norm(R_transpose_R - I, 'fro') > 1e-3  % Frobenius norm per la tolleranza
        disp('La matrice non è ortonormale (R^T * R != Id, tolleranza applicata).');
        disp(R_transpose_R)
        disp('Risultato: Non è una matrice di rotazione.');
        return;
    else
        disp('La matrice è ortonormale.');
    end

    % Check if the determinant is 1
    det_R = det(R);
    if abs(det_R - 1) > 1e-3  % tolleranza per approssimazioni numeriche
        disp(['Il determinante della matrice è ', num2str(det_R), ', non è uguale a 1.']);
        disp('Risultato: Non è una matrice di rotazione.');
        return;
    end
    
    % If all checks pass, the matrix is a rotation matrix
    disp('Risultato: È una matrice di rotazione.');

end


%--------------------------------------------------------------------------------------------------------------------------


%%In the section below, enter the code to solve the exercise.

syms q_1 q_2 l1 l2
N = 2
DHTABLE =  [pi -l1 0 q_1; 
           -pi/2 l2 0 q_2 ];

[T0N, p, n, s, a] = direct_kinematics(N, DHTABLE)
pretty(T0N)

