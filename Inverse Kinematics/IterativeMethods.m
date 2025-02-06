clear variables;
clc;

% Gradient Descent function
function [x, iter] = gradient_descent(f, y, x0, delta, alpha, tol, max_iter, verbose)
    fprintf('\n------------GRADIENT DESCENT METHOD------------\n');
    x = x0;
    for iter = 1:max_iter
        error = f(x) - y;
        norm_error = norm(error);
        J = numJ(f, x, delta);  % Compute Jacobian
        
        % Compute singularity measures
        det_J = det(J);
        [~, S, V] = svd(J);
        min_singular_value = min(diag(S));
        null_space_basis = V(:, end);  % Last column of V is the null space basis
        
        % Check if error is in the null space of J
        is_error_in_null_space = norm(J * null_space_basis) < tol && norm(null_space_basis' * error) < tol;
        
        if verbose
            q_str = sprintf('%f, ', x);  % Create a dynamic string for q
            q_str = q_str(1:end-2);  % Remove trailing comma and space
            fprintf('Iteration %d: q = [%s], norm(error) = %f\n', iter, q_str, norm_error);
            fprintf('  Determinant of J = %f, Min singular value of J = %f\n', det_J, min_singular_value);
            %old
            % fprintf('Iteration %d: q = [%f, %f], norm(error) = %f\n', iter, x(1), x(2), norm_error);
            % fprintf('  Determinant of J = %f, Min singular value of J = %f\n', det_J, min_singular_value);
            if is_error_in_null_space
                fprintf('  Error vector is in the null space of the Jacobian.\n');
            end
        end
        
        % Terminate if error is in the null space and at a singularity
        if is_error_in_null_space && min_singular_value < tol
            warning('Gradient Descent stopped: error vector is in the null space of the Jacobian at a singular configuration.');
            return;
        end
        
        % Terminate if error is within tolerance
        if norm_error < tol
            return;
        end
        
        % Gradient Descent update
        x = x - alpha * J' * error;
    end
    warning('Gradient Descent did not converge within the maximum iterations');
end

% Newton's Method function
function [x, iter] = newton_method(f, y, x0, delta, tol, max_iter, verbose)
    fprintf('\n------------NEWTON METHOD------------\n');
    x = x0;
    for iter = 1:max_iter
        error = f(x) - y;
        norm_error = norm(error);
        J = numJ(f, x, delta);  % Compute Jacobian
        
        % Compute singularity measures
        det_J = det(J);
        [~, S, ~] = svd(J);
        min_singular_value = min(diag(S));
        
        if verbose
            q_str = sprintf('%f, ', x);  % Create a dynamic string for q
            q_str = q_str(1:end-2);  % Remove trailing comma and space
            fprintf('Iteration %d: q = [%s], norm(error) = %.3e\n', iter, q_str, norm_error);
            fprintf('  Determinant of J = %f, Min singular value of J = %f\n', det_J, min_singular_value);
        end       
        if norm_error < tol
            if verbose
                fprintf('Newton Method converged at iteration %d\n', iter);
            end
            return;
        end
        x = x - J \ error;  % Newton's Method update using Jacobian inverse
    end
    disp("!!!!!!!!!!!!!!!")
    disp("-----METODO NEWTON NON è RIUSCITO A CONVERGERE ENTRO IL NUMERO MASSIMO DI ITERAZIONI----")
    disp("!!!!!!!!!!!!!!!")
    warning('Newton Method did not converge within the maximum iterations');
end

% Numerical Jacobian function
function J = numJ(f, x, delta)
    n = numel(x);
    fx = f(x);
    m = numel(fx);
    J = zeros(m, n);
    for i = 1:n
        x_forward = x;
        x_backward = x;
        x_forward(i) = x_forward(i) + delta;
        x_backward(i) = x_backward(i) - delta;
        J(:, i) = (f(x_forward) - f(x_backward)) / (2 * delta);  % Central difference
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the vector-valued function
% l1 = 0.5;
% l2 = 0.4;

% direct kinamatics
f = @(q) [0.5*cos(q(1))+0.5*cos(q(1) + q(2))*cos(q(3));
          0.5*sin(q(1))+0.5*sin(q(1) + q(2))*cos(q(3));
          0.5 + 0.5*sin(q(3))
         ];
% Target value for which we want to find the inverse
% P DESIDERED
y = [0.3; -0.3; 0.7];

% Initial guess
% q0 initial configuration
x0 = [-pi/4; pi/4; pi/4];

% Parameters
delta = 1e-4;  % Small step size for numerical Jacobian
alpha = 0.05;  % Learning rate
tol = 1e-3;    % EPS = Tolerance for convergence -> NORM OF THE ERROR CARTESIAN ---> norm(P_DESIDERED - P_CURRENT)
max_iter = 10;  % Maximum number of iterations
verbose = true;    % Enable detailed output

% Call Gradient Descent
[x_gd, iter_gd] = gradient_descent(f, y, x0, delta, alpha, tol, max_iter, verbose);
% Display Gradient Descent result
fprintf('Gradient Descent converged in %d iterations\n', iter_gd);
y_str = sprintf('%f; ', y);  % Convert 'y' into a semicolon-separated string
y_str = y_str(1:end-2);      % Remove the trailing semicolon and space
x_gd_str = sprintf('%f; ', x_gd);  % Convert 'x_gd' into a semicolon-separated string
x_gd_str = x_gd_str(1:end-2);      % Remove the trailing semicolon and space
fprintf('The inverse solution of f(x) = [%s] is approximately x_sol = [%s]\n', y_str, x_gd_str);
y_GD_sol = f(x_gd);
y_GD_sol_str = sprintf('%f; ', y_GD_sol);  % Convert 'y_GD_sol' into a semicolon-separated string
y_GD_sol_str = y_GD_sol_str(1:end-2);      % Remove the trailing semicolon and space
fprintf('f(x_sol) = [%s]\n', y_GD_sol_str);


% Call Newton's Method
[x_nm, iter_nm] = newton_method(f, y, x0, delta, tol, max_iter, verbose);
% Display Newton's Method result
fprintf('Newton Method HAS STOPPED in %d iterations\n', iter_nm);
y_str = sprintf('%f; ', y);  % Convert 'y' into a semicolon-separated string
y_str = y_str(1:end-2);      % Remove the trailing semicolon and space
x_nm_str = sprintf('%f; ', x_nm);  % Convert 'x_nm' into a semicolon-separated string
x_nm_str = x_nm_str(1:end-2);      % Remove the trailing semicolon and space
fprintf('The inverse solution of f(x) = [%s] is approximately x_sol = [%s]\n', y_str, x_nm_str);
y_NM_sol = f(x_nm);
y_NM_sol_str = sprintf('%f; ', y_NM_sol);  % Convert 'y_NM_sol' into a semicolon-separated string
y_NM_sol_str = y_NM_sol_str(1:end-2);      % Remove the trailing semicolon and space
fprintf('f(x_sol) = [%s]\n', y_NM_sol_str);