function [percent_matrix_qr, percent_matrix_lu] = test_abs_error_random(n,densities,rconds)
    % TEST_ABS_ERROR_RANDOM for a fixed dimension n,density percentage,
    % and condition number, 100 random matrices are generated and the
    % percentage error for successful calculator (absolute value
    % approximation) by taking the eigenvalue decomposition output of the
    % sign to be the 'exact' value the percent_correct_matrix rows
    % represent the various density values while the columns represent the
    % varying reciprocal condition numbers. The outputs percent_matrix_qr
    % and percent_matrix_lu hold as entries in each row and column a
    % percentage that represents how many of the 100 runs succeeded in an
    % absolute error of 0.
    %   [percent_matrix_qr, percent_matrix_lu] = test_abs_error_random(n,densities,rconds)
    numb_densities = length(densities);
    numb_rconds = length(rconds);
    percent_matrix_qr = zeros(numb_rconds+1,numb_densities+1);
    percent_matrix_qr(2:numb_rconds+1,1) = rconds';
    percent_matrix_lu = zeros(numb_rconds+1,numb_densities+1);
    percent_matrix_lu(2:numb_rconds+1,1) = rconds';
    
    number_of_runs = 100;
    for d_ind = 1:length(densities)
        for rc_ind = 1:length(rconds)
            total_correct_qr = 0;
            total_correct_lu = 0;
            R_densities = zeros(number_of_runs);
            parfor k = 1:number_of_runs
                % Generate random matrix
                U = sprand(n,n,densities(d_ind),1);

                % We need to create a block diagonal with some 1x1 and 2x2
                % blocks.

                % Generate values to go in 1x1 blocks and enforce condition
                % number by forcing two values to be 1 and rcond and every
                % other 1x1 block being a value rcond<= j <= 1.
                list_ev = rand(1,n/2);
                list_ev(1) = rconds(rc_ind);
                list_ev(n/2) = 1;

                % generates a list to apply random negatives to the
                % positive numbers in list_ev
                rand_sign = randi([1,2],1,n/2);
                rand_sign = 2*rand_sign - 3*spones(n/2);
                list_ev = rand_sign.*list_ev;
                % Compute the actual sign based on construction and
                sign_actual = prod(rand_sign);
                % Generates 1x1 block diagonal
                D_1 = diag(list_ev);
                % Constants that will determine the 2x2 blocks
                list_bev = linspace(rcond(rc_ind),1,n/4);
                % Generates blcok matrix with 2x2 blocks
                D_2 = sparse(n/2,n/2);
                for l = 1:n/4
                    D_2(2*l-1,2*l) = list_bev(l);
                    D_2(2*l, 2*l-1) = -list_bev(l);
                end
                % Generates the complete block diagonal which includes 1x1
                % and 2x2 blocks
                D = blkdiag(D_1,D_2);
                R = U*D*ctranspose(U);
                R_densities(k) = nnz(R)/(n*n);
                
                % Compute sign(det(A))
                sign_qr = sign_det(R,'qr');
                sign_lu = sign_det(R,'lu');

                % Store value for later use in plot
                abs_error_qr = abs(sign_qr - sign_actual);
                abs_error_lu = abs(sign_lu - sign_actual);

                % Update number of correct absolute differences
                if abs_error_qr == 0
                    total_correct_qr = total_correct_qr + 1;
                end
                if abs_error_lu == 0
                    total_correct_lu = total_correct_lu + 1;
                end
            end
            average_density = sum(R_densities,'all')/length(R_densities);
            percent_matrix_qr(1,d_ind+1) = average_density;
            percent_matrix_lu(1,d_ind+1) = average_density;
            percent_matrix_qr(rc_ind+1,d_ind+1) = total_correct_qr/number_of_runs;
            percent_matrix_lu(rc_ind+1,d_ind+1) = total_correct_lu/number_of_runs;
        end
    end

    writematrix(percent_matrix_lu,'percent_matrix_lu_'+string(n)+'.csv')
    writematrix(percent_matrix_qr,'percent_matrix_qr_'+string(n)+'.csv')
end