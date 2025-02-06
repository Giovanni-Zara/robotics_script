function isFeasible = checkFeasibility(rangeMatrix, testVector)
    %CHECKFEASIBILITY Verifies if a vector is feasible based on a range matrix.
    % INPUTS:
    % rangeMatrix - A matrix representing the range.
    % testVector - The vector to be checked.
    % OUTPUT:
    % isFeasible - A logical value (true if feasible, false otherwise).
    
    % Check dimensions
    [rows, ~] = size(rangeMatrix);
    if length(testVector) ~= rows
        error('---ATTENZIONE--- The length of testVector must match the number of rows in rangeMatrix.');
    end
    
    % Append the test vector as the last column of rangeMatrix
    augmentedMatrix = [rangeMatrix, testVector];
    
    % Compute the rank of the original matrix
    originalRank = rank(rangeMatrix);
    
    % Compute the rank of the augmented matrix
    augmentedRank = rank(augmentedMatrix);
    
    % Check if the rank has changed
    isFeasible = (originalRank == augmentedRank);

    if isFeasible == 1
        disp("Yes it is feasible");
    else
        disp("NO it isn't feasible");
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%syms L q3;
syms gam;
rangeMatrix = [ 1                   0;
                0                   1;
                cos(gam)/sin(gam)   0];

testVector = [0;0;1];

is_feasible = checkFeasibility(rangeMatrix, testVector);
