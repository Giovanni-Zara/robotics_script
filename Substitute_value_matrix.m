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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% syms q1 q2 q3 L q1_dot q2_dot q3_dot;
% J = [- q3_dot*cos(q1 + q3) - q1_dot*(cos(q1 + q3) + q2*cos(q1)) - q2_dot*sin(q1), -q1_dot*sin(q1), -cos(q1 + q3)*(q1_dot + q3_dot);
%   q2_dot*cos(q1) - q1_dot*(sin(q1 + q3) + q2*sin(q1)) - q3_dot*sin(q1 + q3),  q1_dot*cos(q1), -sin(q1 + q3)*(q1_dot + q3_dot);
%                                                                           0,               0,                               0];
% values = struct('q1', pi/2, 'q2', 1, 'q3', 0, 'q1_dot', 1, 'q2_dot', -1, 'q3_dot', -1, 'L', 1);
% substitutedMatrix = substituteMatrix(J, values);
% disp('Substituted Matrix:');
% disp(double(substitutedMatrix));

clear variables
syms q1 q2 l1 l2 real;

% ee pos and orientation -> direct kinematics
J = [- l2*sin(q1 + q2) - l1*sin(q1), -l2*sin(q1 + q2);
  l2*cos(q1 + q2) + l1*cos(q1),  l2*cos(q1 + q2)];
  % End-effector position

values = struct('q1', 0, 'q2', 0, 'l1', 2, 'l2', 1);
substitutedMatrix = substituteMatrix(J, values);
disp('Substituted Matrix:');
disp(double(substitutedMatrix));