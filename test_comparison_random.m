function [percent_matrix_qr, percent_matrix_lu, percent_matrix_qrlu] = test_comparison_random(n,densities,rconds)
    % TEST_COMPARISON_RANDOM
    % for a fixed dimension n,density percentage, and condition number,
    % 100 random matrices are generated and the approximate percentage
    % absolute error for successful calculator (absolute value
    % approximation) by taking the eigenvalue decomposition output of the
    % sign to be the 'exact' value the percent_correct_matrix rows
    % represent the various density values while the columns represent the
    % varying reciprocal condition numbers or more
    % specifically the percent of a successful computation.
    %   [percent_matrix_qr, percent_matrix_lu, percent_matrix_qrlu] = test_comparison_random(n,densities,rconds)
    % Note: This doesn't give a good indication of absolute errors as
    % reciprocal condition number close to zero also generates errors in an
    % eigenvalue decomposition and thus the comparison between the result
    % from using eigs vs lu or qr is not correct.
    numb_densities = length(densities);
    numb_rconds = length(rconds);
    percent_matrix_qr = zeros(numb_rconds+1,numb_densities+1);
    percent_matrix_qr(2:numb_rconds+1,1) = rconds';
    percent_matrix_qr(1,2:numb_densities+1) = densities;
    percent_matrix_lu = zeros(numb_rconds+1,numb_densities+1);
    percent_matrix_lu(2:numb_rconds+1,1) = rconds';
    percent_matrix_lu(1,2:numb_densities+1) = densities;
    percent_matrix_qrlu = zeros(numb_rconds+1,numb_densities+1);
    percent_matrix_qrlu(2:numb_rconds+1,1) = rconds';
    percent_matrix_qrlu(1,2:numb_densities+1) = densities;

    for d_ind = 1:length(densities)
        for rc_ind = 1:length(rconds)
            total_correct_qr = 0;
            total_correct_lu = 0;
            total_correct_qrlu = 0;
            for k = 1:100
                % Generate random matrix
                R = sprand(n,n,densities(d_ind),rconds(rc_ind));
                
                % Compute sign(det(A))
                sign_eigs = sign_det(R,'eigs');
                sign_qr = sign_det(R,'qr');
                sign_lu = sign_det(R,'lu');

                % Store value for later use in plot
                abs_error_qr = abs(sign_eigs - sign_qr);
                abs_error_lu = abs(sign_eigs - sign_lu);
                abs_error_qrlu = abs(sign_qr - sign_lu);

                % Update number of correct absolute differences
                if abs_error_qr == 0
                    total_correct_qr = total_correct_qr + 1;
                end
                if abs_error_lu == 0
                    total_correct_lu = total_correct_lu + 1;
                end
                if abs_error_qrlu == 0
                    total_correct_qrlu = total_correct_qrlu + 1;
                end
            percent_matrix_qr(rc_ind+1,d_ind+1) = total_correct_qr/100;
            percent_matrix_lu(rc_ind+1,d_ind+1) = total_correct_lu/100;
            percent_matrix_qrlu(rc_ind+1,d_ind+1) = total_correct_qrlu/100;
            end
        end
    end
    writematrix(percent_matrix_lu,'percent_matrix_lu_'+string(n)+'.csv')
    writematrix(percent_matrix_qr,'percent_matrix_qr_'+string(n)+'.csv')
    writematrix(percent_matrix_qrlu,'percent_matrix_qrlu_'+string(n)+'.csv')
end