syms c a b w t;

eq_1 = 2*a*w*cos(2*t*w)*cos(t);
der_1 = diff(eq_1, t);
eq_2 = b*w*cos(t*w)*sin(t);
der_2 = diff(eq_2, t);

disp([der_1; der_2] );