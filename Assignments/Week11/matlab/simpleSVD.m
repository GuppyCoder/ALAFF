function [U, S, V] = simpleSVD(B)
    % B is an mxn matrix, and this function computes its SVD
    
    % Step 1: Calculate B'B
    BtB = B' * B;
    
    % Step 2: Calculate eigenvalues of B'B
    eigenvalues = eig(BtB);
    
    % Step 3: The singular values are the square roots of the eigenvalues
    singularValues = sqrt(eigenvalues);
    
    % Sort singular values in descending order
    singularValues = sort(singularValues, 'descend');
    
    % Step 4: Calculate U, S, V matrices
    S = diag(singularValues);
    
    % V are the eigenvectors of B'B
    [V, ~] = eig(BtB);
    
    % U is calculated as B * V * inv(S)
    U = B * V / S;
    
    % Step 5: Display results
    disp('Matrix B:');
    disp(B);
    disp('Singular values matrix S:');
    disp(S);
    disp('Orthogonal matrix U:');
    disp(U);
    disp('Orthogonal matrix V:');
    disp(V);
    
    % Reconstruct matrix
    reconstructedB = U * S * V';
    disp('Reconstructed B from U, S, V:');
    disp(reconstructedB);
end
