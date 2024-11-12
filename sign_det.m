function sgn = sign_det(A, method)
    % SIGN_DET computes the sign of the determinant of a sparse matrix A
    % and method.
    %
    %   sgn = sign_det(A) computes the sign of the determinant of a sparse
    %   matrix A via LU decomposition
    %
    %   sgn = sign_det(A, method) computes the sign of the determinant of
    %   a sparse matrix A with the provided method:
    %   "lu" (LU Decomposition), "qr" (QR Decomposition), or "eigs"
    %   (Eigenvalue decomposition).
    
    % 'Quick check' Shifted Greshgorin Circle Method for symmetric matrices
    % can be used. If successful compute sign determinant and exit
    if method == "gs"
        [success, sgn] = sign_greshgorin(A);
    % If unsuccesful continue with factorization method.
    elseif method == "eigs"
        [r,~] = size(A);
        % Compute Sparse Eigenvalues
        eig_vals = eigs(A,r,'smallestreal');

        % Compute sign of det
        sgn = sign_by_eig_vector(eig_vals);

    elseif method == "lu"
        % Perform Sparse LU Decomposition
        [L, U, P, Q] = lu(A, 'vector');

        % Compute negative entries on L and U
        sign_L = sign_triangular(L);
        sign_U = sign_triangular(U);

        % Compute sing of permutation matrices
        sign_P = sign_perm_vector(P);
        sign_Q = sign_perm_vector(Q);

        % Compute sign of det
        sgn = sign_L * sign_U * sign_P * sign_Q;

    elseif method == "qr"
        % Perform Sparse QR Decomposition
        %[Q,R,P] = qr(A,"vector");
        [Q,R,P] = spqr(A,struct('Q','Householder','permutation','vector')) ;
        
        % Compute negative entries of R
        sign_R = sign_triangular(R);

        % Compute det of Q
        %sign_Q = sign(det(Q));
        % This is the number of hoseholder reflections applied in spqr
        [~, numb_hh_reflects] = size(Q.Tau);
        sign_QP = sign_perm_vector(Q.P);
        sign_Q = (-1)^numb_hh_reflects * sign_QP;

        % Compute sign of permutation matrices
        sign_P = sign_perm_vector(P);

        % Compte sign of det
        sgn = sign_R*sign_Q*sign_P;
    end
end