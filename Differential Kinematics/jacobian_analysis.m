function jacobian_analysis(J, variables)
    % Performs the full analysis of a Jacobian, square or not,
    % displaying: 
    % - Determinant
    % - Singularity
    % - Rank
    % - Nullspace
    % - Range
    %
    % Parameters:
    % J = Jacobian matrix
    % variables = list of symbolic variables (only the ones used to differentiate
    % the Jacobian, don't put constants in there)

    [m, n] = size(J);
    if m ~= n % If m not equal to n (joints) | If m == n, then it's square.
        if n - m >= 2
            % Compute the determinant
            determinant_Jacobian = simplify(det(J.' * J));
        
            disp('Determinant of the Jacobian (J.T * J) --> Since its not square:');
            disp(simplify(determinant_Jacobian));
        
            disp("You might want to check when the Jacobian loses full rank.")
            disp("For doing so, you have to solve the Determinant above for the values which make it == 0.")
            disp("Example: solve(simplify(det(J.' * J)) == 0, variables);")
        elseif m < n
            for i = 1:n
                % Create a submatrix of the original non-square matrix by
                % removing one column at a time
                Jcopy = J;
                Jcopy(:, i) = [];
    
                % Determinant of the i-th submatrix
                fprintf('Determinant of the submatrix with column %d removed: \n', i);
                disp(simplify(det(Jcopy)));
    
                % Compute the singularity (equal the determinant to zero)
                fprintf('Singularity condition of the %s -th submatrix: \n', mat2str(i));
                disp(simplify(det(Jcopy) == 0));
            end
        else 
            % Create all submatrices generated from removing any m-n rows
            if m == 6 && n == 4
                % Create the submatrices by creating all the 4x4 submatrices
                k = 0;
                for i = 1:m
                    % Create a submatrix of the original non-square matrix by
                    % removing one row at a time
                    Jcopy = J;
                    Jcopy(i, :) = [];
                    for j = 1:m-1
                        k = k + 1;
                        Jcopy2 = Jcopy;
                        Jcopy2(j, :) = [];
                        % Determinant of the i-th submatrix
                        fprintf('Determinant of the submatrix with row %d and %d removed: \n', i, j);
                        disp(simplify(det(Jcopy2)));
                    end
                end
            else
                disp('The Jacobian you have parsed has more rows than columns and is not a geometric Jacobian.')
                % Compute the determinant
                determinant_Jacobian = simplify(det(J.' * J));
        
                disp('Determinant of the Jacobian:');
                disp(simplify(determinant_Jacobian));
        
                disp("You might want to check when the Jacobian loses full rank.")
                disp("For doing so, you have to solve the Determinant above for the values which make it == 0.")
                disp("Example: solve(simplify(det(J.' * J)) == 0, variables);")
            end
        end
    else
        % Compute the determinant
        determinant_Jacobian = simplify(det(J));

        disp('Determinant of the Jacobian:');
        disp(simplify(determinant_Jacobian));

        n_var = length(variables);
        % Solve for all possible values of the symbols that make the determinant zero
        if n_var == 2
            [q1_sol, q2_sol, ~, ~] = solve(determinant_Jacobian == 0, variables, 'ReturnConditions', true);
            solutions = [q1_sol, q2_sol];
        elseif n_var == 3
            [q1_sol, q2_sol, q3_sol, ~, ~] = solve(determinant_Jacobian == 0, variables, 'ReturnConditions', true);
            solutions = [q1_sol, q2_sol, q3_sol];
        elseif n_var == 4
            [q1_sol, q2_sol, q3_sol, q4_sol, ~, ~] = solve(determinant_Jacobian == 0, variables, 'ReturnConditions', true);
            solutions = [q1_sol, q2_sol, q3_sol, q4_sol];
        else
            % Handle cases where n_var is not 2, 3, or 4
            solutions = solve(determinant_Jacobian == 0, variables, 'ReturnConditions', true);
        end

        % Display solutions only if they exist
        if exist('solutions', 'var') && ~isempty(solutions)
            disp('All possible solutions for the symbols that make the determinant zero:');
            display(solutions);
        else
            disp('No solutions found for the symbols that make the determinant zero.');
        end
    end

    % Compute the rank of the Jacobian
    rank_J = rank(J);
    
    % Display the rank of the Jacobian
    disp('Rank of the Jacobian:');
    disp(rank_J);
    
    % Compute the nullspace of the Jacobian
    nullspace_J = null(J);
    
    % Display the nullspace of the Jacobian
    disp('Nullspace of the Jacobian:');
    disp(simplify(nullspace_J));

    % Compute the range of the Jacobian
    range = simplify(orth(sym(J)));

    % Display the range
    disp('Range of the Jacobian:');
    disp(range);

    disp('I will display the basis of the range now. -- Warning! This result might be wrong, please check by hand if possible.')
    % Display the basis of the range
    [~, m] = size(range);
    for i = 1:m
        [numerator, ~] = numden(range(:, i));
        fprintf('%d basis\n', i);
        display(simplify(numerator, 'Steps', 50));
    end
    
    disp("--------------------------COMPLEMENTARY (J TRANSPOSED)-------------------------")
    % Compute the complementary nullspace of the Jacobian
    nullspace_J_complementary = null(J.');
    if isempty(nullspace_J_complementary)
        disp('The complementary nullspace is empty.');
    else
        % Display the complementary nullspace of the Jacobian
        disp('Nullspace of the Jacobian (complementary):');
        disp(simplify(nullspace_J_complementary));
    end
    disp("-------------------------------------------------------------------------------")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

N = 3;
syms q1 q2 q3 a1 a2 a3 d1;
J = [  0,  sin(q1), sin(q1);
         0, -cos(q1), -cos(q1);
        1,   0, 0];
 

variables = [q1,q2 q3];
jacobian_analysis(J, variables);

