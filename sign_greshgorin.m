function sgn = sign_greshgorin(A)
    % SIGN_GRESHGORIN Attempts to compute the sign of the determinant via
    % greshgorin circle, however this function only works if the input
    % matrix A is diagonal dominant. It will throw an error if the matrix
    % is not diagonal dominant.
    success_flag = true;
    eig_val_counter = 0;
    [r,~] = size(A);
    % run through diagonals and count number of
    for k = 1:r
        % If diagonal entry is negative.
        if A(k,k)<0
            % Check if upper bound of GS is positive
            radius = sum(abs(A(:,k))) - abs(A(k,k));
            if A(k,k)+radius >= 0
                success_flag = false;
                break
            end
            eig_val_counter = eig_val_counter + 1;
        % If diagonal entry is zero
        elseif A(k,k) == 0
            success_flag = false;

        % If diagonal entry is positive
        else
            % Check if lower bound of GS is negative
            radius = sum(abs(A(:,k))) - abs(A(k,k));
            if A(k,k)-radius <=0
                success_flag = false;
            end
        end
    end
    if success_flag == true
        sgn = (-1)^eig_val_counter;
    else
        error('Matrix is numerically diagonal dominant')
    end
end