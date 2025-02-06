%% Compute the geometric Jacobian
function J = compute_geometric_jacobian(p_vec, z_vec, joints_str)
    JP = [];
    JO = [];
    N = length(joints_str);

    for i = 1:N
        p_i = p_vec(:, i);
        z_i = z_vec(:, i);
        if joints_str(i) == 'R'
            JP = [JP, cross(z_i, p_vec(:, end) - p_i)];
            JO = [JO, z_i];
        else
            JP = [JP, z_i];
            JO = [JO, [0; 0; 0]];
        end
    end

    J = [JP; JO];
    if isa(J, 'sym')
        J = simplify(J);
    end
    
    %disp("Geometric Jacobian matrix:");
    %disp(J);
end

%% Compute the geometric Jacobian for a transformation matrix
function J = transformation_geometric_jacobian(J, R, r)
    M1 = [R, zeros(3); zeros(3), R];
    M2 = [eye(3), -skew(r); zeros(3), eye(3)];
    J = M1 * M2 * J;
end

%% Analyze singularities based on Jacobian determinant
function det_J = analyze_singularities(J, joints_var)
    
    %check if the jacobian is square
    [rows, cols] = size(J);
    JJ = 0;

    if rows == cols
        if isa(JJ, 'sym')
            JJ = simplify(J);
        end
        
    else 
        if rows < cols
            if isa(J * J', 'sym')
                JJ = simplify(J * J');
            end
            fprintf("Jacobian matrix is not square (%d x %d), using J * J'\n", rows, cols);
        else
            if isa(J' * J, 'sym')
                JJ = simplify(J' * J);
            end
            fprintf("Jacobian matrix is not square (%d x %d), using J' * J\n", rows, cols);
        end
    end
    
    if isa(det(JJ), 'sym')
        det_J = simplify(det(JJ));
    else
        det_J = det(JJ);
    end


    %disp("Determinant of the Jacobian:");
    %disp(det_J);

    singular_points = solve(det_J == 0, joints_var);

    if isempty(singular_points)
        disp("No singular configurations found.");
    else 
        disp("Singular configurations: ");
        disp(singular_points);
    end

end

%% Compute the inverse of a transformation matrix
function inverseT = invT(T)
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    inverseT = [R' -R'*p; 0 0 0 1];
end

 %% Compute the skew-symmetric matrix of a vector
function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

%% Build transformation matrices for each link
function A = build_transformation_matrices(DHTABLE)
    N = size(DHTABLE, 1);
    A = cell(1, N);
    for i = 1:N
        alpha = DHTABLE(i, 1);
        a = DHTABLE(i, 2);
        d = DHTABLE(i, 3);
        theta = DHTABLE(i, 4);

        TDH = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
                sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
                  0             sin(alpha)             cos(alpha)            d;
                  0               0                      0                   1];

        A{i} = TDH;
    end
end

%% Perform direct kinematics
function [p_vec, z_vec, T0N] = direct_kinematics(A)
    T = eye(4);
    N = length(A);
    p_i = T(1:3, 4);
    z_i = T(1:3, 3);
    disp("p_0 = [" + join(string(p_i), "; ") + "];");
    disp("z_0 = [" + join(string(z_i), "; ") + "];");
    p_vec = [p_i];
    z_vec = [z_i];

    for i = 1:N
        T = T * A{i};
        if isa(T, 'sym')
            T = simplify(T);
        end
        % disp p_i and z_i
        p_i = T(1:3, 4);
        z_i = T(1:3, 3);
        disp("p_" + i + " = [" + join(string(p_i), "; ") + "];");
        disp("z_" + i + " = [" + join(string(z_i), "; ") + "];");
        p_vec = [p_vec, p_i];
        z_vec = [z_vec, z_i];
    end
    T0N = T;
    fprintf("\n\n__________________________________\n\n");
    disp("Final Transformation Matrix T0N:");
    disp(T0N);
end

function analyze_robot(DH_table, joint_types, joint_vars)
    % Generalized function to analyze robot kinematics
    % Inputs:
    %   - DH_table: Denavit-Hartenberg table (Nx4 matrix)
    %   - joint_types: String indicating joint types ('R' or 'P') (e.g., 'RRR')
    %   - joint_vars: Symbolic variables for joint angles/positions

    % Input validation
    N = size(DH_table, 1); % Number of joints
    assert(N == length(joint_types), "Mismatch between DH table and joint types.");
    assert(N == length(joint_vars), "Mismatch between DH table and joint variables.");

    % Build transformation matrices
    A = build_transformation_matrices(DH_table);

    % Perform direct kinematics
    [p_vec, z_vec, T0N] = direct_kinematics(A);

    % Compute Jacobian
    J_geom = compute_geometric_jacobian(p_vec, z_vec, joint_types);

    % Analyze singularities
    %% da rivedereeeee
    if N < 7 
        detJ = analyze_singularities(J_geom(1:N, :), joint_vars);
    else 
        detJ = analyze_singularities(J_geom(1:6, :), joint_vars);
    end 

    % Display results
    disp("Geometric Jacobian:");
    if isa(J_geom, 'sym')
        disp(simplify(J_geom));
    else
        disp(J_geom);
    end

    disp("Determinant of the Jacobian:");
    if isa(detJ, 'sym')
        disp(simplify(detJ));
    else
        disp(detJ);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear variables;

N = 3;
syms q1 q2 q3 a1 d1 a2 a3;
DHTABLE =  [ pi/2 a1 d1 q1;
                0 a2  0 q2;
             pi/2 a3  0 q3
            ];


joint_types = 'RRR';
joint_vars = [q1, q2, q3];
analyze_robot(DHTABLE, joint_types, joint_vars);

