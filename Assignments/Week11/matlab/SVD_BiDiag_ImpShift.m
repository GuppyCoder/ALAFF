function [S, U, V] = SVD_BiDiag_ImpShift(B)
    [m, n] = size(B);
    U = eye(m);
    V = eye(n);

    % Step 1: Initial Bidiagonalization using Givens rotations
    [B_bidiagonal, U] = applyGivens(B);
    
    % Extract diagonal and superdiagonal elements
    d = diag(B_bidiagonal);
    e = diag(B_bidiagonal, 1);

    % Tolerance and maximum iterations for convergence
    tol = 1e-10;
    maxIter = 1000;
    iter = 0;

    % Step 2: Apply Wilkinson's shift iteratively
    while max(abs(e)) > tol && iter < maxIter
        % Calculate Wilkinson's shift
        mu = wilkinsonShift(d, e);

        % Apply the shift to the bottom 2x2 submatrix
        d(end-1:end) = d(end-1:end) - mu;

        % Reapply Givens rotations to maintain bidiagonal structure
        [B_bidiagonal, U] = applyGivens(B_bidiagonal);

        % Update the diagonal and superdiagonal elements
        d = diag(B_bidiagonal);
        e = diag(B_bidiagonal, 1);

        iter = iter + 1;
    end

    % Step 3: Final SVD on the refined bidiagonal matrix
    [U_bidiagonal, S_diag, V_bidiagonal] = simpleSVD(B_bidiagonal);

    % Step 4: Update U and V
    U = U * U_bidiagonal;
    V = V * V_bidiagonal;

    % Step 5: Extract singular values as a vector
    S = diag(S_diag);
end
