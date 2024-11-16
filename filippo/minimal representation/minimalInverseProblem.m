function res = minimalInverseProblem(R, sequence)
    res = {};
    tol = 1e-6;  % Definisci una tolleranza per il confronto


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

    switch sequence
        case "XYX"
            % cos(theta) in posizione (1, 1)
            % (1, 2), (1, 3)
            ctheta = R(1, 1);
            stheta = sqrt( (R(1, 2))^2 + (R(1, 3))^2 );
            theta1 = atan2(stheta, ctheta);
            theta2 = atan2(-stheta, ctheta);
            
            if stheta > tol    % stheta ~= 0
                % phi
                sphi1 = R(2, 1) / stheta;
                sphi2 = - sphi1;
                cphi1 = - R(3, 1) / stheta;
                cphi2 = - cphi1;
                phi1 = atan2(sphi1, cphi1);
                phi2 = atan2(sphi2, cphi2);
                % psi
                spsi1 = R(1, 2) / stheta;
                spsi2 = - spsi1;
                cpsi1 = - R(1, 3) / stheta;
                cpsi2 = - cpsi1;
                psi1 = atan2(spsi1, cpsi1);
                psi2 = atan2(spsi2, cpsi2);
                
                res{1} = {"ok", phi1, theta1, psi1};
                res{2} = {"ok", phi2, theta2, psi2};
            else
                if theta1 < tol || theta2 < tol    % theta == 0
                    c13 = R(2, 2);
                    s13 = - R(2, 3);
                    sum13 = atan2(s13, c13);
                    res{1} = {"singular0", theta1, sum13};
                    res{2} = {"singular0", theta2, sum13};
                elseif theta1 - pi < tol || theta2 - pi < tol    % theta == pi
                    c13 = R(2, 2);
                    s13 = R(2, 3);
                    diff13 = atan2(s13, c13);
                    res{1} = {"singularpi", theta1, diff13};
                    res{2} = {"singularpi", theta2, diff13};
                end
            end
        case "XYZ"
            % sin(theta) in posizione (1, 3)
        case "XZX"
            % cos(theta) in posizione (1, 1)
        case "XZY"
            % -sin(theta) in posizione (1, 2)
        case "YXY"
            
        case "YXZ"
            
        case "YZX"
            
        case "YZY"
            
        case "ZXY"
            
        case "ZXZ"
            
        case "ZYX"
            
        case "ZYZ"
            
        otherwise

    end

    
end