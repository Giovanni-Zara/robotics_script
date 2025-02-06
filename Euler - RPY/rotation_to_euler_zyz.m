function to_return = rotation_to_euler_zyz(R)
    % rotation_to_euler: Calcola gli angoli di Euler ZYZ dalla matrice di rotazione R.
    % 
    % Inputs:
    %   R - Matrice di rotazione 3x3
    % 
    % Outputs:
    %   uno, due, tre - Angoli di Euler calcolati dalla matrice di rotazione

    %%THE FOLLOWING PART IS FIXED, NO NEED TO EDIT
    % Define expressions for trigonometric components
    c2 = R(3,3);
    s2 = sqrt(R(1,3)^2 + R(2,3)^2);
    
    % Check if s2 is non-zero
    if s2 ~= 0
        disp("-----CASO REGOLARE ESISTONO DUE SOLUZIONI-----");
        disp("-+-+-+-+-+-SOLUZIONE POSITIVA-+-+-+-+-+-")
        c3 = -R(3,1) / s2;
        s3 = R(3,2) / s2;
        c1 = R(1,3) / s2;
        s1 = R(2,3) / s2;
        
        % Display trigonometric components (without pretty for now)
        disp('Trigonometric components based on matrix R:  (s2!=0) ');
        disp(['c2 = r33 = ', num2str(c2)]);
        %c2 = str2double(c2) * pi / 180
        disp(['s2 = sqrt( (r13)^2 + (r23)^2 = )', num2str(s2)]);
        disp(['c3 = -r31/s2 = ', num2str(c3)]);
        disp(['s3 = r32/s2 = ', num2str(s3)]);
        disp(['c1 = r13/s2 = ', num2str(c1)]);
        disp(['s1 = r23/s2 = ', num2str(s1)]);

        uno = atan2(s1,c1);
        due = atan2(s2,c2);
        tre = atan2(s3,c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);

        disp("-+-+-+-+-+-SOLUZIONE NEGATIVA-+-+-+-+-+-")
        s2 = -s2;
        c3 = -R(3,1) / s2;
        s3 = R(3,2) / s2;
        c1 = R(1,3) / s2;
        s1 = R(2,3) / s2;
        
        % Display trigonometric components (without pretty for now)
        disp('Trigonometric components based on matrix R:  (s2!=0) ');
        disp(['c2 = r33 = ', num2str(c2)]);
        %c2 = str2double(c2) * pi / 180
        disp(['s2 = sqrt( (r13)^2 + (r23)^2 = )', num2str(s2)]);
        disp(['c3 = -r31/s2 = ', num2str(c3)]);
        disp(['s3 = r32/s2 = ', num2str(s3)]);
        disp(['c1 = r13/s2 = ', num2str(c1)]);
        disp(['s1 = r23/s2 = ', num2str(s1)]);

        uno = atan2(s1,c1);
        due = atan2(s2,c2);
        tre = atan2(s3,c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);
        
    % Handle special cases when s2 is zero
    elseif s2 == 0 && c2 > 0
        disp('Special case: angolo = 0');
        disp("--------conosci solo sen/cos della somma di angoli--------");
        disp(['cos(uno + tre) = r11 = ', num2str(R(1,1))]);
        disp(['sin(uno + tre) = r21 = ', num2str(R(2,1))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = 0; % Il caso speciale implica angolo zero
        tre = NaN;
        
    elseif s2 == 0 && c2 < 0
        disp('Special case: angolo = pi/-pi');
        disp("--------conosci solo sen/cos della diff di angoli--------");
        disp(['-cos(uno - tre) = r11 = ', num2str(R(1,1))]);
        disp(['-sin(uno - tre) = r21 = ', num2str(R(2,1))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = pi; % Il caso speciale implica angolo zero
        tre = NaN;
    end
    
    to_return = [uno, due, tre];
end

