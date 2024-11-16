function res = axisangleInverseProblem(R)
    res = {};

    % VERIFICA SE R E' UNA ROTATION MATRIX

    tol = 1e-6;  % Definisci una tolleranza per il confronto
    c1 = R(:, 1);
    c2 = R(:, 2);
    c3 = R(:, 3);
    out_prod_condition_1 = abs(c1'*c2)<tol;
    out_prod_condition_2 = abs(c1'*c3)<tol;
    out_prod_condition_3 = abs(c2'*c3)<tol;
    %fprintf("ort12: %d\tort13: %d\tort23: %d\n", out_prod_condition_1, out_prod_condition_2, out_prod_condition_3);
    norm_condition = abs(norm(c1)-1)<tol && abs(norm(c2)-1)<tol && abs(norm(c3)-1)<tol;
    %fprintf("norm1: %d\tnorm2: %d\tnorm3: %d\nnorm_condition: %d\n", norm(c1), norm(c2), norm(c3), norm_condition);
    det_condition = abs(det(R)-1)<tol;
    %fprintf("determinant: %d\tdetcondition: %d\n", det(R), det_condition);

    if ~(out_prod_condition_1 && out_prod_condition_2 && out_prod_condition_3 && norm_condition && det_condition)
        fprintf("La matrice in input NON è una rotation matrix\n");
        res{1}=-1;
        return;
    end
    

    % CALCOLO DELLE SOLUZIONI

    x = ( (R(1, 2)-R(2, 1))^2 + (R(1, 3)-R(3, 1))^2 + (R(2, 3)-R(3, 2))^2 );
    if x < tol   % sin(theta) = 0, caso singolare

        if abs(trace(R)-3) < tol    % theta = 0
            % nessuna soluzione
            fprintf("Non esiste soluzione per r (rotation axis undefined)\n");
            res{1}=-1;
            return;

        else    % theta = +- pi
            % Prima soluzione
            theta1 = pi;
            r1 = 0;
            rx_vals = [sqrt((R(1,1) + 1) / 2), -sqrt((R(1,1) + 1) / 2)];
            ry_vals = [sqrt((R(2,2) + 1) / 2), -sqrt((R(2,2) + 1) / 2)];
            rz_vals = [sqrt((R(3,3) + 1) / 2), -sqrt((R(3,3) + 1) / 2)];

            found = false;
            for rx = rx_vals
                if found, break; end
                for ry = ry_vals
                    if found, break; end
                    for rz = rz_vals
                        cond1 = abs(rx*ry - R(1,2)/2) < tol;
                        cond2 = abs(rx*rz - R(1,3)/2) < tol;
                        cond3 = abs(ry*rz - R(2,3)/2) < tol;
                        if cond1 && cond2 && cond3
                            r1 = [rx; ry; rz];
                            found = true;
                            break
                        end
                    end
                end
            end
            if ~found
                fprintf("Errore, caso singolare, non trovata soluzione\n");
                return;
            end
            res{1} = {theta1, r1};

            % Seconda soluzione
            theta2 = -pi;
            r2 = -r1;
            res{2} = {theta2, r2};
        end

    else
        % Prima soluzione
        theta1 = atan2(sqrt(x), trace(R)-1);
        r1 = 1 / (2 * sin(theta1)) * [R(3, 2) - R(2, 3); R(1, 3) - R(3, 1); R(2, 1) - R(1, 2)];
        res{1} = {theta1, r1};
        
        % Seconda soluzione
        theta2 = -theta1;
        r2 = -r1;
        res{2} = {theta2, r2};
    end

end
