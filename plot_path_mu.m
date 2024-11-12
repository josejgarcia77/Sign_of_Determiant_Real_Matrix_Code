function plot_path_mu(num_of_points)
    % PLOT_PATH_MU plots the topological invariant with the three methods
    % eigs, lu, and qr of a kitaev chain as we
    % change the parameter mu.
    inputs = linspace(-4,4,num_of_points);
    outputs_eigs = zeros(num_of_points,1);
    outputs_lu = zeros(num_of_points,1);
    outputs_qr = zeros(num_of_points,1);

    for index = 1:num_of_points
        mu = inputs(index);
        [X,H_mu] = kitaev_chain(25,2,2,mu,1);
        A = X+ 1i*H_mu;
        outputs_eigs(index) = sign_det(A,'eigs');
        outputs_lu(index) = sign_det(A,'lu');
        outputs_qr(index) = sign_det(A,'qr');
    end
    plot(inputs,outputs_qr,'b--o')
    hold on
    disp(outputs_lu)
    plot(inputs,outputs_lu,'k-^')
    hold on
    plot(inputs,outputs_eigs)
    hold off
    title('Topological Invariant')
    xlabel('\mu')
    ylabel('sign(det(X+iH_\mu))')
end