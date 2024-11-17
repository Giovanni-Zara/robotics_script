%% Moltiplicazione matrici rotazione
syms uno due tre 
Rx = [1,  0,         0;
      0,  cos(uno), -sin(uno);
      0,  sin(uno), cos(uno)];
  
Ry = [cos(due),      0,     sin(due);
      0,            1,      0;
      -sin(due),   0, cos(due)];
  
Rz = [cos(tre), -sin(tre), 0;
       sin(tre), cos(tre),  0;
       0,        0,         1];


disp("UNO->ANGOLO DI X; DUE->ANGOLO DI Y; TRE->ANGOLO DI Z")
%% edit here
R = Rz * Ry