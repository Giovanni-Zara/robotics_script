function pseudo_inv = computePseudoInverse(A)
    % Check if the matrix is square
    [m, n] = size(A);
    isSquare = (m == n);
    
    % Check if the matrix is numeric or symbolic
    isNumeric = isnumeric(A);
    isSymbolic = isa(A, 'sym');
    
    % Compute the pseudo-inverse
    if isNumeric
        pseudo_inv = pinv(A);
    elseif isSymbolic
        pseudo_inv = simplify(pinv(A));
    else
        error('Input matrix must be either numeric or symbolic.');
    end
    
    % If the matrix is square and non-singular, suggest using the inverse
    if isSquare
        if isNumeric
            detA = det(A);
            if detA ~= 0
                fprintf('The matrix is square and non-singular. You could also use the inverse.\n');
                disp("THE INVERSE HAS BEEN RETURNED");
                invA = inv(A);
                % Optionally, you can return the inverse instead of the pseudo-inverse
                pseudo_inv = invA;
            end
        elseif isSymbolic
            detA = det(A);
            if detA ~= sym(0)
                fprintf('The matrix is square and non-singular. You could also use the inverse.\n');
                invA = simplify(inv(A));
                disp("THE INVERSE HAS BEEN RETURNED");
                % Optionally, you can return the inverse instead of the pseudo-inverse
                pseudo_inv = invA;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

clear variables
syms a1 a3 a4 real;
J = [    0     0;
     3     1];
 
J_pinv = computePseudoInverse(J);
disp(J_pinv*[0;-1])
