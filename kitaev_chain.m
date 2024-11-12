function [X,H] = kitaev_chain(N, t, delta, mu, sf_spacing)
    % KITAEV_CHAIN  Generates a kitaev chain real space position operator
    % and hamiltonian in the majorana basis.
    % 
    %   [X,H] = kitaev_chain(N) generates kitaev chain
    %   with nearest neighbor hoping t=2, superconducting pairing term
    %   delta=2, and chemical potential mu=0. The position operator will
    %   be defined with unit spacing between the spinless fermions in 1d.
    %
    %   [X,H] = kitaev_chain(N, t, delta, mu) generates kitaev chain
    %   with nearest neighbor hoping t, superconducting pairing term
    %   delta, and chemical potential mu. The position operator will
    %   be defined with unit spacing between the spinless fermions in 1d.
    %   
    %   [X,H] = kitaev_chain(N, t, delta, mu, sf_spacing) generates kitaev,
    %   chain with nearest neighbor hoping t, superconducting pairing term
    %   delta, and chemical potential mu. The position operator will
    %   be defined with sf_spacing constant between the spinless fermions
    %   in 1d.
    %   
    % For more information, see <a href="kitaevchain: 
    % web('https://crangi.github.io/post/kitaev_chain/')">
    % A basic introduction to Majorana Edge modes in a Kitaev Chain</a>.
    
    % Initiate Sparse Matrices for Position and Hamiltonian
    X = sparse(2*N,2*N);
    H = sparse(2*N,2*N);

    switch nargin
        case 4
            % Default position spacing.
            sf_spacing = 1;
        case 1
            % Default parameters and spacing.
            t = 2;
            delta = 2;
            mu = 0;
            sf_spacing = 1;
    end

    % Initiate constants to be used repeatedly.
    Jx = 0.5*(t-delta);
    Jy = 0.5*(t+delta);

    % Construct real space Hamiltonian.
    for n = 1:N-1
        H(2*n,2*n+1) = Jx;
        H(2*n+1,2*n) = -Jx;
        H(2*n-1,2*n+2) = -Jy;
        H(2*n+2,2*n-1) = Jy;
        H(2*n-1,2*n) = mu;
        H(2*n,2*n-1) = -mu;
    end

    H(2*N-1,2*N) = mu;
    H(2*N,2*N-1) = -mu;

    H = 1i*H;
    
    % Construct position operator with unit*sf_spacing spacing between
    % spinless fermions and centered around 0.

    % Checks to see if we have even or odd number of spinless fermions to
    % keep the lattice centered around zero.
    if mod(N,2) == 0
        start = sf_spacing*(-N/2 + 1/2);
    else
        start = sf_spacing*(-floor(N/2));
    end

    % Creates diagonal matrix with corresponding one dimensional position
    % information.
    for n = 0:N-1
        X(2*n+1,2*n+1) = start + sf_spacing*n;
        X(2*n+2,2*n+2) = start + sf_spacing*n;
end