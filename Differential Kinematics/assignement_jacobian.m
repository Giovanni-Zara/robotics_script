function [rank_J, range_space, null_space, dim_range, dim_null] = analyze_jacobian(J, conditions)
    % analyze_jacobian: Analyzes the Jacobian matrix after substituting conditions
    % Inputs:
    %   J: Jacobian matrix (symbolic or numeric)
    %   conditions: A cell array of conditions to substitute into the Jacobian
    % Outputs:
    %   rank_J: Rank of the Jacobian after substitution
    %   range_space: Basis for the range space of the Jacobian
    %   null_space: Basis for the null space of the Jacobian
    %   dim_range: Dimension of the range space
    %   dim_null: Dimension of the null space

    % Step 1: Substitute the conditions iteratively until no further changes occur
    J_sub = J;
    J_prev = J_sub; % Store the previous state of the Jacobian
    while true
        for i = 1:length(conditions)
            J_sub = subs(J_sub, conditions{i}{1}, conditions{i}{2});
        end
        % Simplify the Jacobian to ensure all substitutions are applied
        J_sub = simplify(J_sub);
        % Check if the Jacobian has stopped changing
        if isequal(J_sub, J_prev)
            break; % Exit the loop if no further changes occur
        end
        J_prev = J_sub; % Update the previous state
    end

    % Step 2: Check if the Jacobian is fully numeric
    if isa(J_sub, 'sym')
        % If the Jacobian still contains symbolic variables, compute rank, range space, and null space symbolically
        rank_J = rank(J_sub);
        range_space = colspace(J_sub); % Use colspace for symbolic range space
        null_space = null(J_sub); % Use null for symbolic null space
    else
        % If the Jacobian is fully numeric, proceed with numeric computations
        rank_J = rank(J_sub);
        range_space = orth(J_sub); % Use orth for numeric range space
        null_space = null(J_sub); % Use null for numeric null space
    end

    % Step 3: Compute the dimensions of the range space and null space
    dim_range = rank_J; % Dimension of the range space is equal to the rank
    dim_null = size(J_sub, 2) - rank_J; % Dimension of the null space is (number of columns - rank)

    % Display the results
    fprintf('Substituted Jacobian:\n');
    disp(J_sub);
    fprintf('Rank of the Jacobian: %d\n', rank_J);
    fprintf('Dimension of the range space: %d\n', dim_range);
    fprintf('Dimension of the null space: %d\n', dim_null);
    fprintf('Range space basis:\n');
    disp(range_space);
    fprintf('Null space basis:\n');
    disp(null_space);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% syms q1 q2 q3
% J = [q1 q2 q3; q2 q3 q1; q3 q1 q2];
% conditions = {{q3, 0}, {q2, -sin(q3)}};
% 
% [rank_J, range_space, null_space] = analyze_jacobian(J, conditions);


syms q1 q2 L1 L2 real;

J = [- L2*sin(q1 + q2) - L1*sin(q1), -L2*sin(q1 + q2);
      L2*cos(q1 + q2) + L1*cos(q1),  L2*cos(q1 + q2)];

conditions = {{q1, pi}, {q2, pi}};

[rank_J, range_space, null_space] = analyze_jacobian(J, conditions);
