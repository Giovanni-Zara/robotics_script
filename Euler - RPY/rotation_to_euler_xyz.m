function to_return = rotation_to_euler_xyz(R)
    % rotation_to_euler: Calcola gli angoli di Euler da una matrice di rotazione R.
    % La combinazione di Euler in questo caso è XYZ.
    % 
    % Inputs:
    %   R - Matrice di rotazione 3x3
    % 
    % Outputs:
    %   uno, due, tre - Angoli di Euler calcolati dalla matrice di rotazione
        
    % Calcola i componenti trigonometrici dalla matrice di rotazione R
    s2 = R(1,3);
    c2 = sqrt(R(2,3)^2 + R(3,3)^2);
    
    disp('UNO-> angolo 1 rotazione, DUE-> angolo 2 rotazione, TRE->angolo 3 rotazione');
    
    % Controlla se c2 è diverso da zero
    if c2 ~= 0
        disp("-----CASO REGOLARE ESISTONO DUE SOLUZIONI-----");
        disp("-+-+-+-+-+-SOLUZIONE POSITIVA-+-+-+-+-+-")
        s3 = -R(1,2) / c2;
        c3 = R(1,1) / c2;
        s1 = -R(2,3) / c2;
        c1 = R(3,3) / c2;

        % Mostra i componenti trigonometrici
        disp('Componenti trigonometriche basate sulla matrice R: (c2 != 0)');
        disp(['c2 = sqrt( (r23)^2 + (r33)^2 ) = ', num2str(c2)]);
        disp(['s2 = r13 = ', num2str(s2)]);
        disp(['c3 = r11 / c2 = ', num2str(c3)]);
        disp(['s3 = -r12 / c2 = ', num2str(s3)]);
        disp(['c1 = r33 / c2 = ', num2str(c1)]);
        disp(['s1 = -r23 / c2 = ', num2str(s1)]);
        
        % Calcola gli angoli di Euler
        uno = atan2(s1, c1);
        due = atan2(s2, c2);
        tre = atan2(s3, c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);

        disp("-+-+-+-+-+-SOLUZIONE NEGATIVA-+-+-+-+-+-")
        c2 = -c2;
        s3 = -R(1,2) / c2;
        c3 = R(1,1) / c2;
        s1 = -R(2,3) / c2;
        c1 = R(3,3) / c2;

        % Mostra i componenti trigonometrici
        disp('Componenti trigonometriche basate sulla matrice R: (c2 != 0)');
        disp(['c2 = sqrt( (r23)^2 + (r33)^2 ) = ', num2str(c2)]);
        disp(['s2 = r13 = ', num2str(s2)]);
        disp(['c3 = r11 / c2 = ', num2str(c3)]);
        disp(['s3 = -r12 / c2 = ', num2str(s3)]);
        disp(['c1 = r33 / c2 = ', num2str(c1)]);
        disp(['s1 = -r23 / c2 = ', num2str(s1)]);
        
        % Calcola gli angoli di Euler
        uno = atan2(s1, c1);
        due = atan2(s2, c2);
        tre = atan2(s3, c3);
        disp(['UNO = ', num2str(uno)]);
        disp(['DUE = ', num2str(due)]);
        disp(['TRE = ', num2str(tre)]);
        
    elseif c2 == 0 && s2 > 0
        disp('Special case: angolo due = +pi/2');
        disp("--------conosci solo sen/cos della somma di angoli--------");
        disp(['cos(uno + tre) = r22 = ', num2str(R(2,2))]);
        disp(['sin(uno + tre) = r21 = ', num2str(R(2,1))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = pi/2; % Il caso speciale implica angolo zero
        tre = NaN;

    elseif c2 == 0 && s2 < 0
        disp('Special case: angolo due = -pi/2');
        disp("--------conosci solo sen/cos della diff di angoli--------");
        disp(['cos(uno - tre) = r22 = ', num2str(R(2,2))]);
        disp(['-sin(uno - tre) = r21 = ', num2str(R(2,1))]);
        uno = NaN; % Indica che il valore non è calcolabile
        due = -pi/2; % Angolo massimo o minimo
        tre = NaN;
    end

    to_return = [uno, due, tre];
end


