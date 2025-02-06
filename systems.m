% Define symbolic variables
syms q_2 px py

% Given equation
eq1 = (2*cos(q_2))/5 + 41/100 == px^2 + py^2;

% Solve for q_2
solution = solve(eq1, q_2);

% Display the solution
disp('The solution for q_2 is:');
disp(vpa(solution, 5)); % Display the solution with 5 decimal places