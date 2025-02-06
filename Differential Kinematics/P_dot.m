function [p_dot, p_ddot] = computeDerivatives(p)
    % computeDerivatives Computes the first and second derivatives of a vector p.
    %
    % Inputs:
    % - p: A symbolic vector representing direct kinematics.
    %
    % Outputs:
    % - p_dot: The first derivative of p with respect to time.
    % - p_ddot: The second derivative of p with respect to time.

    % Ensure the input is symbolic
    if ~isa(p, 'sym')
        error('Input vector p must be symbolic.');
    end

    % Define symbolic time variable
    syms t

    % Compute the first derivative of p with respect to t
    p_dot = diff(p, t);

    % Compute the second derivative of p with respect to t
    p_ddot = diff(p_dot, t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms t q1(t) q2(t);

% Example kinematic vector in 2D or 3D
p = [q2(t)*cos(q1(t));
    q2(t)*sin(q1(t))];

[p_dot, p_ddot] = computeDerivatives(p);

disp('First derivative (velocity):');
disp(p_dot);

disp('Second derivative (acceleration):');
disp(p_ddot);


