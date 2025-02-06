function to_return = rotation_to_euler_xzy(R)
    % rotation_to_euler_xzy: Calcola gli angoli di Euler dalla matrice di rotazione R con configurazione XZY.
    % 
    % Inputs:
    %   R - Matrice di rotazione 3x3
    % 
    % Outputs:
    %   uno, due, tre - Angoli di Euler calcolati dalla matrice di rotazione con configurazione XZY

    %% THE FOLLOWING PART IS FIXED
    % Define expressions for trigonometric components
    s2 = -R(1,2);
    % c2 = sqrt(R(2,2)^2 + R(3,2)^2); % cambiato a seguitoesame luglio 2024
    c2 = sqrt(R(1,1)^2 + R(1,3)^2);

    % Check if c2 is non-zero
    if c2 ~= 0
        disp("-----CASO REGOLARE ESISTONO DUE SOLUZIONI-----");
        disp("-+-+-+-+-+-SOLUZIONE POSITIVA-+-+-+-+-+-")
        % s3 = R(3,1) / c2;
        s3 = R(1, 3) / c2; % cambiato guardando sol
        c3 = R(1,1) / c2;
        s1 = R(3,2) / c2;
        c1 = R(2,2) / c2;
        
        % Display trigonometric components (without pretty for now)
        disp('Trigonometric components based on matrix R:  (s2!=0) ');
        disp(['c2 = sqrt( (r22)^2 + (r32)^2 ) = ', num2str(c2)]);
        disp(['s2 = -r12 = ', num2str(s2)]);
        disp(['c3 = r11/c2 = ', num2str(c3)]);
        disp(['s3 = r31/c2 = ', num2str(s3)]);
        disp(['c1 = r22/c2 = ', num2str(c1)]);
        disp(['s1 = r32/c2 = ', num2str(s1)]);

        uno = atan2(s1, c1);
        due = atan2(s2, c2);
        tre = atan2(s3, c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);

        disp("-+-+-+-+-+-SOLUZIONE NEGATIVA-+-+-+-+-+-")
        c2 = -c2;
        % s3 = R(3,1) / c2;
        s3 = R(1, 3) / c2; % cambiato guardando sol
        c3 = R(1,1) / c2;
        s1 = R(3,2) / c2;
        c1 = R(2,2) / c2;
        
        % Display trigonometric components (without pretty for now)
        disp('Trigonometric components based on matrix R:  (s2!=0) ');
        disp(['c2 = sqrt( (r22)^2 + (r32)^2 ) = ', num2str(c2)]);
        disp(['s2 = -r12 = ', num2str(s2)]);
        disp(['c3 = r11/c2 = ', num2str(c3)]);
        disp(['s3 = r31/c2 = ', num2str(s3)]);
        disp(['c1 = r22/c2 = ', num2str(c1)]);
        disp(['s1 = r32/c2 = ', num2str(s1)]);

        uno = atan2(s1, c1);
        due = atan2(s2, c2);
        tre = atan2(s3, c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);
        
    % Handle special cases when c2 is zero
    elseif c2 == 0 && s2 > 0
        disp('Special case: angolo = pi/2');
        disp("--------conosci solo sen/cos della DIFF di angoli--------");        
        disp(['cos(uno - tre) = r21 = ', num2str(R(2,1))]);
        disp(['-sin(uno - tre) = r23 = ', num2str(R(2,3))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = pi/2; % Il caso speciale implica angolo zero
        tre = NaN;
        
    elseif c2 == 0 && s2 < 0
        disp('Special case: angolo = -pi/2');
        disp("--------conosci solo sen/cos della somma di angoli--------");
        disp(['-sin(uno + tre) = r23 = ', num2str(R(2,3))]);
        disp(['-cos(uno + tre) = r21 = ', num2str(R(2,1))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = -pi/2; % Il caso speciale implica angolo zero
        tre = NaN;
    end
    
    to_return = [uno, due, tre];
end

