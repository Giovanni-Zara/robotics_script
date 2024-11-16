function R = minimalDirectProblem(phi, theta, psi, sequence, fixed_axes)
    % rotation around fixed axes
    if fixed_axes
        % sequenza invertita
        sequence = sequence(end:-1:1);
        % inverte ordine degli angoli
        temp = phi;
        phi = psi;
        psi = temp;
    end

    switch sequence
        case "XYX"
            R1 = rotationX(phi);
            R2 = rotationY(theta);
            R3 = rotationX(psi);
        case "XYZ"
            R1 = rotationX(phi);
            R2 = rotationY(theta);
            R3 = rotationZ(psi);
        case "XZX"
            R1 = rotationX(phi);
            R2 = rotationZ(theta);
            R3 = rotationX(psi);
        case "XZY"
            R1 = rotationX(phi);
            R2 = rotationZ(theta);
            R3 = rotationY(psi);
        case "YXY"
            R1 = rotationY(phi);
            R2 = rotationX(theta);
            R3 = rotationY(psi);
        case "YXZ"
            R1 = rotationY(phi);
            R2 = rotationX(theta);
            R3 = rotationZ(psi);
        case "YZX"
            R1 = rotationY(phi);
            R2 = rotationZ(theta);
            R3 = rotationX(psi);
        case "YZY"
            R1 = rotationY(phi);
            R2 = rotationZ(theta);
            R3 = rotationY(psi);
        case "ZXY"
            R1 = rotationZ(phi);
            R2 = rotationX(theta);
            R3 = rotationY(psi);
        case "ZXZ"
            R1 = rotationZ(phi);
            R2 = rotationX(theta);
            R3 = rotationZ(psi);
        case "ZYX"
            R1 = rotationZ(phi);
            R2 = rotationY(theta);
            R3 = rotationX(psi);
        case "ZYZ"
            R1 = rotationZ(phi);
            R2 = rotationY(theta);
            R3 = rotationZ(psi);
        otherwise
            R1 = -1;
            R2 = -1;
            R3 = -1;
    end
    
    R = R1 * R2 * R3;

end


% funzioni ausiliarie

% rotazione elementare attorno all'asse X di angolo a
function R = rotationX(a)
    R = [1     0        0;
         0     cos(a)   -sin(a);
         0     sin(a)   cos(a)];
end
% rotazione elementare attorno all'asse Y di angolo a
function R = rotationY(a)
    R = [cos(a)     0       sin(a);
         0          1       0;
         -sin(a)    0       cos(a)];
end
% rotazione elementare attorno all'asse Z di angolo a
function R = rotationZ(a)
    R = [cos(a)    -sin(a)     0;
         sin(a)    cos(a)      0;
         0         0           1];
end




