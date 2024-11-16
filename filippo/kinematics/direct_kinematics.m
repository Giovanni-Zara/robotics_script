
clear all
%clc

n=4;    % number of joints

% alfa_i   |   a_i   |   d_i   |   theta_i
DH_matrix = [0      sym('a1')   sym('d1')   sym('q1');
             0      sym('a2')   0           sym('q2');
             0      0           sym('q3')   0;
             pi     0           sym('d4')   sym('q4')];







T = eye(n);
for i = 1:n

    alpha = DH_matrix(i, 1);
    a = DH_matrix(i, 2);
    d = DH_matrix(i, 3);
    theta = DH_matrix(i, 4);

    fprintf("HOM matrix from %d to %d:", i-1, i);
    % homogeneous transformation matrix from RFi-1 to RFi
    A =    [cos(theta)  -sin(theta)*cos(alpha)  sin(theta)*sin(alpha)   a*cos(theta);
            sin(theta)  cos(theta)*cos(alpha)   -cos(theta)*sin(alpha)  a*sin(theta);
            0           sin(alpha)              cos(alpha)              d;
            0           0                       0                       1]

    T = T * A;
    T = simplify(T);    % solo per matrici simboliche

end

fprintf("Homogeneous transformation matrix:\nT =\n\n");
disp(T);

fprintf("\nPosition of RF_n w.r.t. RF_0:");
p = T(1:3,4)












