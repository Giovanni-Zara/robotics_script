%% WITH PHI
function [q1_elbow_up, q2_elbow_up, q1_elbow_down, q2_elbow_down] = computeInverseKinematicsWithOrientation(Px, Py, phi, L1, L2)
    % computeInverseKinematicsWithOrientation - Computes the joint angles q1 and q2 for a planar 2R robot
    % given the end-effector position (Px, Py) and orientation (phi), and link lengths L1 and L2.
    % Supports both numeric and symbolic inputs.
    %
    % Inputs:
    %   Px, Py - Coordinates of the end-effector position (numeric or symbolic)
    %   phi - Orientation of the end-effector (in radians, numeric or symbolic)
    %   L1, L2 - Lengths of the two links (numeric or symbolic)
    %
    % Outputs:
    %   q1_elbow_up, q2_elbow_up - Joint angles for the elbow-up configuration (numeric or symbolic)
    %   q1_elbow_down, q2_elbow_down - Joint angles for the elbow-down configuration (numeric or symbolic)

    fprintf('Starting inverse kinematics computation for Px = %s, Py = %s, phi = %s, L1 = %s, L2 = %s\n', ...
            string(Px), string(Py), string(phi), string(L1), string(L2));

    % Step 1: Compute cos(q2)
    numerator_c2 = Px^2 + Py^2 - L1^2 - L2^2;
    denominator_c2 = 2 * L1 * L2;
    c2 = numerator_c2 / denominator_c2;

    fprintf('Computed cos(q2) = %s\n', string(c2));

    % Check if the point is within the robot's workspace (only for numeric inputs)
    if isnumeric(c2) && abs(c2) > 1
        error('The point (%s, %s) is outside the robot workspace.', string(Px), string(Py));
    end

    % Step 2: Compute sin(q2) for both elbow-up and elbow-down configurations
    s2_elbow_up = sqrt(1 - c2^2);  % Elbow-up: positive sin(q2)
    s2_elbow_down = -sqrt(1 - c2^2);  % Elbow-down: negative sin(q2)

    fprintf('Computed sin(q2) for elbow-up = %s\n', string(s2_elbow_up));
    fprintf('Computed sin(q2) for elbow-down = %s\n', string(s2_elbow_down));

    % Step 3: Compute q2 for both configurations
    q2_elbow_up = atan2(s2_elbow_up, c2);  % Elbow-up
    q2_elbow_down = atan2(s2_elbow_down, c2);  % Elbow-down

    fprintf('Computed q2 for elbow-up = %s radians\n', string(q2_elbow_up));
    fprintf('Computed q2 for elbow-down = %s radians\n', string(q2_elbow_down));

    % Step 4: Compute q1 for both configurations using q1 = phi - q2
    q1_elbow_up = phi - q2_elbow_up;  % Elbow-up
    q1_elbow_down = phi - q2_elbow_down;  % Elbow-down

    fprintf('Computed q1 for elbow-up = %s radians\n', string(q1_elbow_up));
    fprintf('Computed q1 for elbow-down = %s radians\n', string(q1_elbow_down));

    fprintf('Inverse kinematics computation completed.\n');
end

%% WITHOUT PHI
function [q1_elbow_up, q2_elbow_up, q1_elbow_down, q2_elbow_down] = computeInverseKinematics(Px, Py, L1, L2)
    % computeInverseKinematics - Computes the joint angles q1 and q2 for a planar 2R robot
    % given the end-effector position (Px, Py) and link lengths L1 and L2.
    % Supports both numeric and symbolic inputs.
    %
    % Inputs:
    %   Px, Py - Coordinates of the end-effector position (numeric or symbolic)
    %   L1, L2 - Lengths of the two links (numeric or symbolic)
    %
    % Outputs:
    %   q1_elbow_up, q2_elbow_up - Joint angles for the elbow-up configuration (numeric or symbolic)
    %   q1_elbow_down, q2_elbow_down - Joint angles for the elbow-down configuration (numeric or symbolic)

    fprintf('Starting inverse kinematics computation for Px = %s, Py = %s, L1 = %s, L2 = %s\n', ...
            string(Px), string(Py), string(L1), string(L2));

    % Step 1: Compute cos(q2)
    numerator_c2 = Px^2 + Py^2 - L1^2 - L2^2;
    denominator_c2 = 2 * L1 * L2;
    c2 = numerator_c2 / denominator_c2;

    fprintf('Computed cos(q2) = %s\n', string(c2));

    % Check if the point is within the robot's workspace (only for numeric inputs)
    if isnumeric(c2) && abs(c2) > 1
        error('The point (%s, %s) is outside the robot workspace.', string(Px), string(Py));
    end

    % Step 2: Compute sin(q2) for both elbow-up and elbow-down configurations
    s2_elbow_up = sqrt(1 - c2^2);  % Elbow-up: positive sin(q2)
    s2_elbow_down = -sqrt(1 - c2^2);  % Elbow-down: negative sin(q2)

    fprintf('Computed sin(q2) for elbow-up = %s\n', string(s2_elbow_up));
    fprintf('Computed sin(q2) for elbow-down = %s\n', string(s2_elbow_down));

    % Step 3: Compute q2 for both configurations
    q2_elbow_up = atan2(s2_elbow_up, c2);  % Elbow-up
    q2_elbow_down = atan2(s2_elbow_down, c2);  % Elbow-down

    fprintf('Computed q2 for elbow-up = %s radians\n', string(q2_elbow_up));
    fprintf('Computed q2 for elbow-down = %s radians\n', string(q2_elbow_down));

    % Step 4: Compute q1 for both configurations
    % Elbow-up
    numerator_c1_up = Px * (L1 + L2 * c2) + Py * L2 * s2_elbow_up;
    numerator_s1_up = Py * (L1 + L2 * c2) - Px * L2 * s2_elbow_up;
    denominator = L1^2 + L2^2 + 2 * L1 * L2 * c2;

    c1_up = numerator_c1_up / denominator;
    s1_up = numerator_s1_up / denominator;
    q1_elbow_up = atan2(s1_up, c1_up);

    % Elbow-down
    numerator_c1_down = Px * (L1 + L2 * c2) + Py * L2 * s2_elbow_down;
    numerator_s1_down = Py * (L1 + L2 * c2) - Px * L2 * s2_elbow_down;

    c1_down = numerator_c1_down / denominator;
    s1_down = numerator_s1_down / denominator;
    q1_elbow_down = atan2(s1_down, c1_down);

    fprintf('Computed q1 for elbow-up = %s radians\n', string(q1_elbow_up));
    fprintf('Computed q1 for elbow-down = %s radians\n', string(q1_elbow_down));

    fprintf('Inverse kinematics computation completed.\n');
end

%% MAIN FUNCTION 
function [q1, q2] = mainInverseKinematics(Px, Py, varargin)
    % mainInverseKinematics - Main function to compute inverse kinematics for a planar 2R robot.
    % It chooses between two functions based on the input provided.
    % Supports both numeric and symbolic inputs.
    %
    % Inputs:
    %   Px, Py - Coordinates of the end-effector position (numeric or symbolic)
    %   varargin - Optional inputs:
    %       - If varargin contains phi (orientation), it uses computeInverseKinematicsWithOrientation.
    %       - If varargin contains L1 and L2 (link lengths), it uses computeInverseKinematics.
    %
    % Outputs:
    %   q1, q2 - Joint angles in radians (returns both elbow-up and elbow-down solutions if applicable)

    fprintf('Starting main inverse kinematics computation for Px = %s, Py = %s\n', string(Px), string(Py));

    % Check if the optional input includes phi (orientation)
    if nargin == 5
        % Case 1: Input includes phi (orientation)
        phi = varargin{1};
        L1 = varargin{2};
        L2 = varargin{3};

        disp("----------------------------------------------");
        disp("-----ORIENTATION HAS BEEN PASSED AS INPUT-----");
        disp("----------------------------------------------");
        fprintf('Orientation phi = %s provided. Using computeInverseKinematicsWithOrientation.\n', string(phi));
    
        % Call the function for inverse kinematics with orientation
        [q1_elbow_up, q2_elbow_up, q1_elbow_down, q2_elbow_down] = computeInverseKinematicsWithOrientation(Px, Py, phi, L1, L2);
    
        % Return both solutions (elbow-up and elbow-down)
        q1 = [q1_elbow_up; q1_elbow_down];
        q2 = [q2_elbow_up; q2_elbow_down];
    
    elseif nargin == 4
        % Case 2: Input does not include phi (orientation)
        L1 = varargin{1};
        L2 = varargin{2};

        disp("----------------------------------------------");
        disp("--------------WITHOUT ORIENTATION-------------");
        disp("----------------------------------------------");
        fprintf('No orientation provided. Using computeInverseKinematics.\n');
    
        % Call the function for inverse kinematics without orientation
        [q1_elbow_up, q2_elbow_up, q1_elbow_down, q2_elbow_down] = computeInverseKinematics(Px, Py, L1, L2);
    
        % Return both solutions (elbow-up and elbow-down)
        q1 = [q1_elbow_up; q1_elbow_down];
        q2 = [q2_elbow_up; q2_elbow_down];
    
    else
        error('Invalid number of inputs. Expected 4 (Px, Py, L1, L2) or 5 (Px, Py, phi, L1, L2).');
    end

    fprintf('Main inverse kinematics computation completed.\n');
end

%% use case
clear varibales

px = 0.3708;
py = 0.6739;
phi = pi / 4;
l1 = 0.5;
l2 = 0.4;

% non cambiare mai ordine di questi elementi in input, se non hai phi
% rimuovila da input senza cambiare ordine

% YOU GOT PHI
%[q1, q2] = mainInverseKinematics(px, py, phi, l1, l2);

% no PHI
[q1, q2] = mainInverseKinematics(px, py, l1, l2);