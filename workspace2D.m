syms theta1 theta2 theta3 theta4 theta5 theta6 real
syms d1 d2 d3 d4 real

% DH TABLE,  L'ORDINE DEI PARAMETRI è A,ALPHA,D,THETA!!!!!!!!!!!!!!
test_dh = [18.5 0 42 theta1;
           16 pi 0 theta2;
           0 0 d3 0;
           0 0 7 theta4]

% Parameter ranges
theta1_range = arr2Rad(linspace(0,300, 50));
theta2_range = arr2Rad(linspace(0,300, 50));
d3_range = linspace(0,12, 50);
% Note the specification states 540°, but anything past 360° is redundant
theta4_range = arr2Rad(linspace(0,360, 50));
test_map = containers.Map({'theta1', 'theta2', 'd3','theta4'}, ...
    {theta1_range, theta2_range, d3_range, theta4_range}); 
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