%% Analytical Jacobian Analysis

function detJ = analyze_robot_kinematics(joints, links, position)
    % Function to analyze robot kinematics
    % Inputs:
    %   - joints: symbolic array of joint variables (e.g., [q1, q2, q3])
    %   - links: symbolic array of link lengths (e.g., [l1, l2, l3])
    %   - position: symbolic expression of the end-effector position as a column vector

    % Check compatibility between number of joint variables and link lengths
    if numel(joints) < numel(links)
        error('The number of joints cannot be less than the number of link lengths.');
    elseif numel(joints) > numel(links)
        warning('The number of joints exceeds the number of link lengths. Verify the configuration.');
    end

    % Compute the analytical Jacobian
    J = jacobian(position, joints);
    disp('Analytical Jacobian matrix:');
    disp(J);

    % Compute determinant of the Jacobian (or J*J' for redundancy resolution)
    detJ = simplify(det(J * J'));
    disp('Determinant of the Jacobian:');
    disp(detJ);

    % Compute and display singularities
    % Generalized singularity analysis by solving det(J*J') = 0
    singularities = solve(detJ == 0, joints, 'ReturnConditions', true);
    disp('Singularities:');
    disp(singularities.conditions);

    % Compute the null space of the Jacobian
    nullJ = null(J);
    nullJ = simplify(nullJ);
    disp('Null space of the Jacobian:');
    disp(nullJ);

    % Compute the null space of the transpose of the Jacobian
    nullJT = null(J');
    nullJT = simplify(nullJT);
    disp('Null space of the transpose of the Jacobian:');
    disp(nullJT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;

N = 3;
syms q1 q2 q3 l real;

joints = [q1, q2, q3]; % Joint variables
links = [1, 1, 1];  % Link lengths

% ee pos and orientation -> direct kinematics
position = [ q1 + l*cos(q2) + l*cos(q3+q2);
             l*sin(q2) + l*sin(q3+q2);
             q2 + q3
             ];                            

% Call the function
analyze_robot_kinematics(joints, links, position);



