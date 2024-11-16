clear all
syms phi theta psi


% esempio minimalDirectProblem --------------------------------------------
disp("---------------------------------- minimal representation Direct Problem:")

s = "ZYZ";
fixed = false;
R = minimalDirectProblem(phi, theta, psi, s, fixed);

axes = "moving"; axes(fixed) = "fixed";
fprintf("Rotation matrix R for rotation around %s axes %s:\n", axes, s);
disp(R);
if R == -1
    fprintf("Errore: sequenza di rotazione non valida\n");
end


%{
% esempio minimalInverseProblem -------------------------------------------
disp("---------------------------------- minimal representation Inverse Problem:")
R = [-1 0 0; 0 -1/sqrt(2) -1/sqrt(2); 0 -1/sqrt(2) 1/sqrt(2)];
fprintf("Input:\n\trotation matrix:\n");
disp(R);

res = minimalInverseProblem(R, "XYX");
if res{1}{1} == "ok"
    fprintf("Soluzione1:\n\tphi1 = %f, theta1 = %f, psi1 = %f\n\n", res{1}{2}, res{1}{2}, res{1}{3});
    fprintf("Soluzione2:\n\tphi2 = %f, theta2 = %f, psi2 = %f\n", res{2}{2}, res{2}{2}, res{2}{3});
elseif res{1}{1} == "singular0"
    fprintf("Soluzione:\n\ttheta1 = %f, theta2 = %f, phi+psi = %f\n\n", res{1}{2}, res{2}{2}, res{1}{3});
elseif res{1}{1} == "singularpi"
    fprintf("Soluzione:\n\ttheta1 = %f, theta2 = %f, phi-psi = %f\n\n", res{1}{2}, res{2}{2}, res{1}{3});
else
    disp("Errore");
end
%}



















