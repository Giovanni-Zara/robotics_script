% Define angles as symbolic variables
syms uno due tre angolo
disp('UNO->ANGOLO X, DUE->ANGOLO Y, TRE->ANGOLO Z')

% Define rotation matrices
Rx = [1,  0,         0;
      0,  cos(uno), -sin(uno);
      0,  sin(uno), cos(uno)];
  
Rz = [cos(due), -sin(due), 0;
      sin(due), cos(due),  0;
      0,        0,         1];
  
Rx1 = [1,  0,         0;
      0,  cos(tre), -sin(tre);
      0,  sin(tre), cos(tre)];

% Combined rotation matrix
xzx = Rx * Rz * Rx1;

%% EDIT HERE THE NUMBERS
% Assign numeric values to angles if u need numbers - DEPENDS ON EXERCISE
a = pi/3;  %edit here
b = pi/3;  %edit here
c = pi/3;  %edit here
%%
% Define R matrix using numeric values - IT DEPENDS ON THE EXERCISE U HAVE!
%%For the rotation that u have in the exercises-a, b, c here represents the
%%angles given in the exercise. a->x, b->y, c->z
Rxa = [1,  0,         0;
      0,  cos(a), -sin(a);
      0,  sin(a), cos(a)];
Ryb = [cos(b),      0,     sin(b);
      0,            1,      0;
      -sin(b),       0,    cos(b)];
Rzc =[cos(c), -sin(c), 0;
      sin(c), cos(c),  0;
      0,        0,     1];

%% EDIT HERE FOR rotation matrix that DEPENDS ON THE EXERCISE
R = Rzc * Ryb;      %edit here
%%

%%FOLLOWING PART IS FIXED
disp('matrice Rxzx')
print_matrix(xzx, 2)

s2 = sqrt(R(2,1)^2 + R(3,1)^2);
c2 = R(11);

% Check if c2 is non-zero
if s2 ~= 0
    s3 = R(1,3) / s2;
    c3 = -R(1,2) / s2;
    s1 = R(3,1) / s2;
    c1 = R(2,1) / s2;
    
    % Display trigonometric components (without pretty for now)
    disp('Trigonometric components based on matrix R:  (s2!=0) ');

    disp(['c2 = r11 = ', num2str(c2)]);
    %c2 = str2double(c2) * pi / 180
    disp(['s2 = sqrt( (r21)^2 + (r31)^2 ) = ', num2str(s2)]);

    disp(['c3 = -r12/s2 = ', num2str(c3)]);
    disp(['s3 = -r13/s2 = ', num2str(s3)]);
    disp(['c1 = r21/s2 = ', num2str(c1)]);
    disp(['s1 = r31/s2 = ', num2str(s1)]);
    
% Handle special cases when s2 is zero
elseif s2 == 0 && c2 > 0
    R(2,2) = cos(uno + tre);
    R(3,2) = sin(uno + tre);
    disp('Special case: angolo = 0');
    disp(['r22 = ', num2str(R(2,2))]);
    disp(['r32 = ', num2str(R(3,2))]);
    
elseif s2 == 0 && c2 < 0
    R(2,2) = -cos(uno - tre);
    R(3,2) = -sin(uno - tre);
    disp('Special case: angolo = +/- pi');
    disp(['r22 = ', num2str(R(2,2))]);
    disp(['r32 = ', num2str(R(3,2))]);
end

disp('UNO->ANGOLO X, DUE->ANGOLO Y, TRE->ANGOLO Z')
uno = atan2(s1,c1)
due = atan2(c2,c2)
tre = atan2(s3,c3)



function print_matrix(matrix, precision)
    % print_matrix: Prints a symbolic matrix in a clear format
    %
    % Inputs:
    % - matrix: symbolic matrix to print
    % - precision (optional): number of decimal places to use (if numeric output desired)
    %
    % If precision is not specified, the matrix will be printed symbolically.

    % Check if precision is provided for numeric output
    if nargin < 2
        precision = []; % Default: no numeric approximation
    end
    
    % Display message for the matrix
    disp('Matrix output:')
    
    % Use symbolic pretty printing if precision is not specified
    if isempty(precision)
        pretty(matrix)
    else
        % Use vpa for decimal approximation with specified precision
        matrix_numeric = vpa(matrix, precision);
        disp(matrix_numeric);
    end
end
