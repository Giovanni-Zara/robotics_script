% Define angles as symbolic variables
syms uno due tre angolo 
disp('UNO->ANGOLO X, DUE->ANGOLO Y, TRE->ANGOLO Z')


% Define rotation matrices
Rx = [1,  0,         0;
      0,  cos(uno), -sin(uno);
      0,  sin(uno), cos(uno)];
  
Ry = [cos(due),      0,     sin(due);
      0,            1,      0;
      -sin(due),   0, cos(due)];
  
Rz = [cos(tre), -sin(tre), 0;
       sin(tre), cos(tre),  0;
       0,        0,         1];

% Combined rotation matrix
xyz = Rx * Ry * Rz;

%% EDIT HERE THE NUMBERS
% Assign numeric values to angles if u need numbers - DEPENDS ON EXERCISE
a = 0;
b = -pi/3;
c = pi/4;
%%
% Define R matrix using numeric values - IT DEPENDS ON THE EXERCISE U HAVE!
%%For the rotation that u have in the exercises
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
R = Rzc * Ryb;
%%

%%FOLLOWING PART IS FIXED
disp('matrice Rxyz')
print_matrix(xyz, 2)

s2 = R(1,3);
c2 = sqrt(R(2,3)^2 + R(3,3)^2);

% Check if c2 is non-zero
if c2 ~= 0
    s3 = -R(1,2) / c2;
    c3 = R(1,1) / c2;
    s1 = -R(2,3) / c2;
    c1 = R(3,3) / c2;
    

    % Display trigonometric components (without pretty for now)
    disp('Trigonometric components based on matrix R:  (c2!=0) ');

    disp(['c2 = sqrt( (r23)^2 + (r33)^2 ) = ', num2str(c2)]);
    %c2 = str2double(c2) * pi / 180
    disp(['s2 = r13 = ', num2str(s2)]);

    disp(['c3 = r11/c2 = ', num2str(c3)]);
    disp(['s3 = -r12/c2 = ', num2str(s3)]);
    disp(['c1 = r32/c2 = ', num2str(c1)]);
    disp(['s1 = -r23/c2 = ', num2str(s1)]);
    
% Handle special cases when s2 is zero
elseif c2 == 0 && s2 > 0
    R(2,2) = cos(uno + tre);
    R(2,1) = sin(due + tre);
    disp('Special case: angolo = pi/2');
    disp(['r22 = ', num2str(R(2,2))]);
    disp(['r21 = ', num2str(R(2,1))]);
    
elseif c2 == 0 && s2 < 0
    R(2,2) = cos(uno - tre);
    R(2,1) = sin(uno - tre);
    disp('Special case: angolo = -pi/2');
    disp(['r22 = ', num2str(R(2,2))]);
    disp(['r21 = ', num2str(R(2,1))]);
end

disp('UNO->ANGOLO X, DUE->ANGOLO Y, TRE->ANGOLO Z')
uno = atan2(s1,c1)
due = atan2(s2,c2)
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
