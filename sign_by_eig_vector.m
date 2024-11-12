function sgn = sign_by_eig_vector(eigenvalue_vector)
    % SIGN_BY_EIG_VECTOR calculates the sign of a matrix given an input
    % vector of eigenvalues sorted from smallest to largest.
    % 
    sgn = 1;
    sign_counter = 0;
    % Check if eigenvalue_vector is a column vector, if it is then
    % converts it to a row vector, otherwise no action is done.
    if iscolumn(eigenvalue_vector)
        eigenvalue_vector = eigenvalue_vector';
    end

    % Loop through eigenvalues to determine compute sign break if zero is
    % encountered or if we encounter positive values
    for val = eigenvalue_vector
        if val < 0
            sign_counter = sign_counter + 1;
        elseif val == 0
            sgn = 0;
            break;
        else
            break;
        end
    end
    % If above loops completes and sign is not zero then we compute the
    % sign based on number of negative eigenvalues.
    if not(sgn == 0)
        sgn = (-1)^(sign_counter);
    end
end