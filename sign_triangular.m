function sgn = sign_triangular(T)
    % SIGN_TRIANGULAR computes the sign of a triangular (upper or lower)
    % by counting the number of negative entries on the diagonal.
    sgn = 1;
    [r,~] = size(T);
    sign_counter = 0;
    for k = 1:r
        if T(k,k) < 0
            sign_counter = sign_counter + 1;
        elseif T(k,k) == 0
            sgn = 0;
            break
        end
    end
    if not(sgn==0)
        sgn = (-1)^(sign_counter);
    end
end