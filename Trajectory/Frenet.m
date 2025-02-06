function [t, n, b] = Frenet_axis(p, derive_wrt)
    % computeDerivatives Computes the first and second derivatives of a vector p.
    %
    % Inputs:
    % - p: A symbolic vector representing direct kinematics.
    % - derive_wrt: The variable with respect to which to differentiate.
    %
    % Outputs:
    % - t: The tangent vector.
    % - n: The normal vector.
    % - b: The binormal vector.

    % Ensure the input is symbolic
    if ~isa(p, 'sym')
        error('Input vector p must be symbolic.');
    end
  
    % Compute the first derivative of p with respect to derive_wrt
    p_dot = diff(p, derive_wrt);
    
    % Compute the tangent vector t
    norm_p_dot = norm(p_dot);
    if isa(p_dot, 'sym')
         p_dot = simplify(p_dot);
         norm_p_dot = simplify(norm(p_dot));
    end
    t = p_dot / norm_p_dot;
    disp("the T vector, expressed 1/norm * vect")
    disp(norm_p_dot);
    disp(p_dot);

    disp("---------------------------------------------")
    % Compute the second derivative of p with respect to derive_wrt
    p_ddot = diff(p_dot, derive_wrt);
    disp("p_ddot")
    disp(p_ddot)
    
    % Compute the normal vector n
    t_dot = diff(t, derive_wrt);
    norm_t_dot = norm(t_dot);
    if isa(t_dot, 'sym')
         t_dot = simplify(t_dot);
         norm_t_dot = simplify(norm(t_dot));
    end
    n = t_dot / norm_t_dot;
    disp("the N vector, expressed 1/norm * vect")
    disp(norm_t_dot)
    disp(t_dot)
    
    disp("---------------------------------------------")
    % Compute the binormal vector b
    b = cross(t, n);
    if isa(b, 'sym')
        b = simplify(b);
    end
    disp("the B vector");
    disp(b);
    disp("---------------------------------------------")
    
    % Simplify the vectors if they are symbolic
    if isa(t, 'sym')
        t = simplify(t);
        n = simplify(n);
        b = simplify(b);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define symbolic variables
syms r h s real;

% Example kinematic vector in 3D
p = [r*sin(s);
     h*s;
     r*(1 + cos(s))];

% Define the variable with respect to which to differentiate
derive_wrt = s;

% Compute the Frenet frame
[t, n, b] = Frenet_axis(p, derive_wrt);

% Display the results
disp('Tangent vector t:');
disp(t);
disp('Normal vector n:');
disp(n);
disp('Binormal vector b:');
disp(b);