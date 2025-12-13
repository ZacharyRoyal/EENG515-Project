function coeffs = alt_disc_fourier(samples, metric_def)
    N = size(samples, 2);
    Phi = construct_discrete_sample_matrix(N, metric_def);
    b = Phi * samples';

    % we now have our vector b in the normal equations, we now need to get
    % G, the grammian matrix, invert it, and multiply it by b
    G = Phi' * Phi/N; % we already have phi, calling construct grammian matrix would recalculate it 
    % and waste a bunch of time

    [Q, R] = qr(G);
    coeffs = (inv(R) * Q' * b)';
    % fuck you matlab, inv is better than your wacky division syntax
end