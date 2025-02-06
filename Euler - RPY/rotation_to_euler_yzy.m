function to_return = rotation_to_euler_yzy(R)
    % rotation_to_euler_yzy: Calcola gli angoli di Euler da una matrice di rotazione R.
    % La combinazione di Euler in questo caso è YZY.
    %
    % Inputs:
    % R - Matrice di rotazione 3x3
    %
    % Outputs:
    % uno, due, tre - Angoli di Euler calcolati dalla matrice di rotazione
    
    s2 = sqrt(R(1,2)^2 + R(3,2)^2);
    c2 = R(2,2);
    
    % Check if c2 is non-zero
    if s2 ~= 0
        disp("-----CASO REGOLARE ESISTONO DUE SOLUZIONI-----");
        disp("-+-+-+-+-+-SOLUZIONE POSITIVA-+-+-+-+-+-")
        s3 = R(2,3) / s2;
        c3 = R(2,1) / s2;
        s1 = R(3,2) / s2;
        c1 = -R(1,2) / s2;
        
        % Display trigonometric components (without pretty for now)
        disp('Trigonometric components based on matrix R: (s2!=0) ');
        disp(['c2 = r22 = ', num2str(c2)]);
        disp(['s2 = sqrt( (r12)^2 + (r32)^2 ) = ', num2str(s2)]);
        disp(['c3 = r21/s2 = ', num2str(c3)]);
        disp(['s3 = r23/s2 = ', num2str(s3)]);
        disp(['c1 = -r12/s2 = ', num2str(c1)]);
        disp(['s1 = r32/s2 = ', num2str(s1)]);
    
        disp('UNO-> angolo 1 rotazione, DUE-> angolo 2 rotazione, TRE->angolo 3 rotazione');
        uno = atan2(s1,c1);
        due = atan2(s2,c2);
        tre = atan2(s3,c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);

        disp("-+-+-+-+-+-SOLUZIONE NEGATIVA-+-+-+-+-+-")
        s2 = -s2;
        s3 = R(2,3) / s2;
        c3 = R(2,1) / s2;
        s1 = R(3,2) / s2;
        c1 = -R(1,2) / s2;
        
        % Display trigonometric components (without pretty for now)
        disp('Trigonometric components based on matrix R: (s2!=0) ');
        disp(['c2 = r22 = ', num2str(c2)]);
        disp(['s2 = sqrt( (r12)^2 + (r32)^2 ) = ', num2str(s2)]);
        disp(['c3 = r21/s2 = ', num2str(c3)]);
        disp(['s3 = r23/s2 = ', num2str(s3)]);
        disp(['c1 = -r12/s2 = ', num2str(c1)]);
        disp(['s1 = r32/s2 = ', num2str(s1)]);
    
        disp('UNO-> angolo 1 rotazione, DUE-> angolo 2 rotazione, TRE->angolo 3 rotazione');
        uno = atan2(s1,c1);
        due = atan2(s2,c2);
        tre = atan2(s3,c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);

    % Handle special cases when s2 is zero
    elseif s2 == 0 && c2 > 0
        disp('Special case: angolo due = 0');
        disp("--------conosci solo sen/cos della somma di angoli--------");        
        disp(['cos(uno + tre) = r11 = ', num2str(R(1,1))]);
        disp(['sin(uno + tre) = r13 = ', num2str(R(1,3))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = 0; % Il caso speciale implica angolo zero
        tre = NaN;

    elseif s2 == 0 && c2 < 0
        disp('Special case: angolo due = +/- pi');
        disp("--------conosci solo sen/cos della somma di angoli--------");
        disp(['-cos(uno - tre) = r11 = ', num2str(R(1,1))]);
        disp(['sin(uno - tre) = r13 = ', num2str(R(1,3))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = +pi; % Il caso speciale implica angolo zero
        tre = NaN;
    end
    
    to_return = [uno, due, tre];
end