function R = axisangleDirectProblem(theta, r)
    if norm(r) ~= 1
        % normalizza l'asse r per ottenere un vettore unitario
        r = r / norm(r);
    end
    
    % componenti dell'asse
    rx = r(1);
    ry = r(2);
    rz = r(3);

    c = cos(theta);
    s = sin(theta);

    % rotation matrix
    R = [ rx^2*(1-c)+c,     rx*ry*(1-c)-rz*s,   rx*rz*(1-c)+ry*s;
          rx*ry*(1-c)+rz*s, ry^2*(1-c)+c,       ry*rz*(1-c)-rx*s;
          rx*rz*(1-c)-ry*s, ry*rz*(1-c)+rx*s,   rz^2*(1-c)+c; ];
end
