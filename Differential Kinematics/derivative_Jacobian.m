function substitutedMatrix = substituteMatrix(matrix, values)
    % substituteMatrix Substitutes values into a symbolic matrix.
    %
    % Inputs:
    % - matrix: A symbolic matrix.
    % - values: A structure containing parameter names as fields and values to substitute.
    %
    % Outputs:
    % - substitutedMatrix: The matrix after substituting values.
    
    % Ensure the input matrix is symbolic
    if ~isa(matrix, 'sym')
        error('Input matrix must be symbolic.');
    end
    
    % Convert the structure to symbolic variables and values
    parameters = fieldnames(values); % Get field names
    parameterValues = struct2cell(values); % Get corresponding values
    
    % Substitute the values into the matrix
    substitutedMatrix = subs(matrix, parameters, parameterValues);
end


function simplifiedMatrix = simplifyMatrix(matrix)
    % simplifyMatrix Simplifies a symbolic matrix.
    %
    % Inputs:
    % - matrix: A symbolic matrix to simplify.
    %
    % Outputs:
    % - simplifiedMatrix: The simplified version of the input matrix.
    
    % Ensure the input is symbolic
    if ~isa(matrix, 'sym')
        error('Input matrix must be symbolic.');
    end

    % Apply simplification to each element
    simplifiedMatrix = simplify(matrix);
end



function dJ_dt = compute_jacobian_derivative(J, q, q_dot)
    % Computes the derivative of the Jacobian matrix with respect to joint configurations
    %
    % Inputs:
    %   J     - Symbolic Jacobian matrix (m x n)
    %   q     - Symbolic vector of joint variables (1 x n or n x 1)
    %   q_dot - Symbolic or numeric vector of joint velocities (1 x n or n x 1)
    %
    % Output:
    %   dJ_dt - Derivative of the Jacobian matrix (m x n)
    
    % Initialize the derivative matrix
    [m, n] = size(J);  % Size of the Jacobian matrix
    dJ_dt = sym(zeros(m, n));  % Symbolic matrix of the same size as J

    % Compute the derivative term for each joint
    for j = 1:n
        % Compute the partial derivative of J with respect to q(j)
        dJ_dqj = diff(J, q(j));
        
        % Multiply by the corresponding joint velocity q_dot(j)
        dJ_dt = dJ_dt + dJ_dqj * q_dot(j);
    end
end



% Define symbolic variables
% syms q1 q2 q3 d1 a2 a3 q1_dot q2_dot q3_dot real;
% q = [q1; q2; q3];       % Joint variables
% q_dot = [q1_dot; q2_dot; q3_dot];  % Joint velocities
% 
% % Example Jacobian (3x3 symbolic matrix)
% J = [-sin(q1)*(a3*cos(q2 + q3) + a2*cos(q2)), -cos(q1)*(a3*sin(q2 + q3) + a2*sin(q2)), -a3*sin(q2 + q3)*cos(q1);
%     cos(q1)*(a3*cos(q2 + q3) + a2*cos(q2)), -sin(q1)*(a3*sin(q2 + q3) + a2*sin(q2)), -a3*sin(q2 + q3)*sin(q1);
%                                      0,            a3*cos(q2 + q3) + a2*cos(q2),          a3*cos(q2 + q3)];
% 
% % Compute the derivative of the Jacobian
% dJ_dt = compute_jacobian_derivative(J, q, q_dot);
% 
% % Display the result
% disp('Derivative of the Jacobian (dJ/dt):');
% disp(dJ_dt);
% 
% simpl = simplifyMatrix(dJ_dt);
% disp(simpl);
% 
% values = struct('q1', pi/2, 'q2', pi/4, 'q3', pi/2, 'q1_dot', 1, 'q2_dot', 2, 'q3_dot', -2, 'd1', 5, 'a2', 4, 'a3', 3);
% 
% substitutedMatrix = substituteMatrix(simpl, values);
% 
% disp('Substituted Matrix:');
% disp(double(substitutedMatrix));


syms q1 q2 q1_dot q2_dot real;
q = [q1; q2];       % Joint variables
q_dot = [q1_dot; q2_dot];  % Joint velocities
L = 1;

J = [0, 1;
     1, 0];

% Compute the derivative of the Jacobian
dJ_dt = compute_jacobian_derivative(J, q, q_dot);

% Display the result
disp('Derivative of the Jacobian (dJ/dt):');
disp(dJ_dt);

simpl = simplifyMatrix(dJ_dt);
disp(simpl);



