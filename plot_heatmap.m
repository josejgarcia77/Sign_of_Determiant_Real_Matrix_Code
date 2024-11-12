function plot_heatmap(num_sites,num_of_points)
    % PLOT_HEATMAP Generates the heatmap of a kitaev chain with num_sites
    % being the number of cites in the kitaev chain, num_of_point being the
    % number of points to compute for x and y axis in plot.
    x_inputs = linspace(-num_sites-10,num_sites+10,num_of_points);
    mu_inputs = linspace(-4,4,num_of_points);
    outputs_eigs = zeros(num_of_points,num_of_points);
    outputs_lu = zeros(num_of_points,num_of_points);
    outputs_qr = zeros(num_of_points,num_of_points);

    for row = 1:num_of_points
        for column = 1:num_of_points
            x = x_inputs(row);
            mu = mu_inputs(column);
            [X,H] = kitaev_chain(num_sites,2,2,mu,1);
            A = (X-x*speye(2*num_sites)) + 1i*H;
            %outputs_eigs((2*25)+1-row,column) = sign_det(A,'eigs');
            outputs_lu((num_of_points+1)-row,column) = sign_det(A,'lu');
            %outputs_qr((2*25)+1-row,column) = sign_det(A,'qr');
        end
    end
    CustomXLabels = cat(2,string(-num_sites-10),strings(1,(num_of_points-3)/2),"0",strings(1,(num_of_points-3)/2),string(num_sites+10));
    CustomYLabels = cat(2,"-4",strings(1,(num_of_points-3)/2),"0",strings(1,(num_of_points-3)/2),"4");
    hm = heatmap(outputs_lu);
    hm.GridVisible = 'off';
    hm.XDisplayLabels = CustomXLabels;
    hm.YDisplayLabels = CustomYLabels;
    %hold off
    title('Topological Invariant')
    xlabel('x')
    ylabel('\mu')
end