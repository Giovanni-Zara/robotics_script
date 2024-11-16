%{
    axisangleDirectProblem:
        Given: theta, r
        Find: rotation matrix for rotation of angle theta around generic axis r

    axisangleInverseProblem:
        Given: rotation matrix R
        Find: theta, r
%}

clear all

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ esempio axisangleDirectProblem
disp("---------------------------------- axis/angle Direct Problem:")

syms theta;     % Angolo di rotazione (in radianti)
%theta = pi/4;
r = [1/sqrt(2), -1/sqrt(2), 0];      % Asse di rotazione

%fprintf("Input:\n\ttheta: %f \n\tasse di rotazione: [%f  %f  %f]\n", theta, r(1), r(2), r(3));

% Calcola la matrice di rotazione
fprintf("\nSoluzione:\nrotation matrix:\n");
R = axisangleDirectProblem(theta, r)

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$














% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ esempio axisangleInverseProblem
disp("---------------------------------- axis/angle Inverse Problem:")

fprintf("Input:\n\trotation matrix:\n");
R = [-1 0 0; 0 -1/sqrt(2) -1/sqrt(2); 0 -1/sqrt(2) 1/sqrt(2)]

% calcola asse di rotazone r e angolo theta
res = axisangleInverseProblem(R);

if ~isequal(res{1}, -1)
    fprintf("\nSoluzione 1:\n");
    fprintf("r:\n");
    disp(res{1}{2});
    fprintf("theta:\n");
    disp(res{1}{1});
    fprintf("\nSoluzione 2:\n");
    fprintf("r:\n");
    disp(res{2}{2});
    fprintf("theta:\n");
    disp(res{2}{1});
end

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$












