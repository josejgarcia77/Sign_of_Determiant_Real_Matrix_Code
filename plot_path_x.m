function plot_path_x(num_of_points)
    % PLOT_PATH_X plots the topological invariant with the three methods
    % using eigenvalues, qr, and lu decompositions as we change the
    % position 'probe' of the system
    inputs = linspace(-30,30,num_of_points);
    outputs_eigs = zeros(num_of_points,1);
    outputs_lu = zeros(num_of_points,1);
    outputs_qr = zeros(num_of_points,1);

    for index = 1:num_of_points
        x = inputs(index);
        [X,H] = kitaev_chain(25,2,2,4,1);
        A = (X-x*speye(2*25)) + 1i*H;
        outputs_eigs(index) = sign_det(A,'eigs');
        outputs_lu(index) = sign_det(A,'lu');
        outputs_qr(index) = sign_det(A,'qr');
    end
    plot(inputs,outputs_qr,'b--o')
    hold on
    plot(inputs,outputs_lu,'k-^')
    hold on
    plot(inputs,outputs_eigs)
    hold off
    title('Topological Invariant')
    xlabel('x')
    ylabel('sign(det((X-x)+iH_\mu))')
end