function [B, U] = applyGivens(B)
    [m, n] = size(B);
    U = eye(m); % Initialize U as an identity matrix

    for j = 1:n
        for i = m:-1:(j+1)
            G = eye(m); % Initialize G as an identity matrix
            [c, s] = givensRotation(B(i-1, j), B(i, j));
            G([i-1, i], [i-1, i]) = [c, s; -s, c]; % Apply rotation to relevant submatrix

            % Apply Givens rotation from the left (on B) and accumulate in U
            B = G' * B;
            U = U * G;
        end
    end
end

