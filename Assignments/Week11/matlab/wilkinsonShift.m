function mu = wilkinsonShift(d, e)
    % d: diagonal elements of the bidiagonal matrix
    % e: superdiagonal elements of the bidiagonal matrix
    
    if isscalar(e)
        % When the matrix is 2x2, directly compute the shift
        delta = (d(1) - d(2)) / 2;
        mu = d(2) - sign(delta) * (e(1)^2) / (abs(delta) + sqrt(delta^2 + e(1)^2));
    else
        % General case for larger matrices
        delta = (d(end-1) - d(end)) / 2;
        mu = d(end) - sign(delta) * (e(end-1)^2) / (abs(delta) + sqrt(delta^2 + e(end-1)^2));
    end
end
