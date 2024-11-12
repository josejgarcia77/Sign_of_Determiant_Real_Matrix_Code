function speed_test_random(density, rcond)
    % SPEED_TEST_RANDOM tests the speed of the qr and lu algorithm for
    % given density and reciprocal condition number rcond.
    % For each matrix dimension n = [10e1,10e2,10e3,2x10e3,3x10e3] compute
    % sign(det(A)) 100 times and print results in box plots.
    x = zeros(100,10);
    n = [100, 300, 600, 1000, 1300];
    t_f = 0;
    t_g = 0;
    for k_n = 1:5
        for k = 1:100
            % Generate random matrix
            R = sprand(n(k_n),n(k_n),density,rcond);

            % Compute sign(det(A))
            f = @() sign_det(R,'qr');
            g = @() sign_det(R,'lu');
            t_f = timeit(f);
            t_g = timeit(g);

            % Store value for later use in plot
            x(k,k_n*2-1) = t_f;
            x(k,k_n*2) = t_g;
        end
    end
    figure
    n_labels = {'qr 100', 'lu 100', 'qr 300', 'lu 300', 'qr 600', 'lu 600', 'qr 1000', 'lu 1000', 'qr 1300', 'lu 1300'};
    boxplot(x,'Labels',n_labels)
    title_label = append('Compare Run for various matrix sizes and methods');
    title(title_label);
    xlabel('Method and size of square matrix')
    ylabel('Time in seconds')
end