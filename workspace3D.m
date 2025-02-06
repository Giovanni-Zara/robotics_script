% Symbolic vars
syms theta1 theta2 theta3 theta4 theta5 theta6 real
syms d1 d2 d3 d4 real

% DH parameters     L'ORDINE DEI PARAMETRI è A,ALPHA,D,THETA!!!!!!!!!!!!!!
test_dh = [0 0    0  theta1;
           0 pi/2 d2 pi/2;
           0 0    d3 0]

% Parameter ranges
theta1_range = arr2Rad(linspace(-149.4, 149.4));
d2_range = linspace(0,50);
d3_range = linspace(10, 30);
% Note the specification states 540°, but anything past 360° is redundant
test_map = containers.Map({'theta1', 'd2', 'd3'}, {theta1_range, d2_range, d3_range});
% Workspace plotting function
plot3dworkspace(test_dh, test_map, @get_alternative_dh_transform)



%FUNZIONI - NON MODIFICARE
function out = arr2Rad(A)
    out = arrayfun(@(angle) deg2rad(angle), A);
end

function T = get_alternative_dh_transform(a,alpha,d,theta)
T = [cos(theta) -cos(alpha)*sin(theta) sin(alpha)*sin(theta) a*cos(theta)
     sin(theta) cos(alpha)*cos(theta) -sin(alpha)*sin(theta) a*sin(theta)
     0 sin(alpha) cos(alpha) d
     0 0 0 1];
end