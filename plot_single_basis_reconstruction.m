function plot_single_basis_reconstruction(samples, metric_def)
    % samples: signal vector
    % metric_def: struct defining the non-L2 metric (e.g. struct('p',6))

    N = numel(samples);
    freqs = 0:(N-1);

    % Build basis matrix for this metric
    Phi = build_alt_basis(N, metric_def);

    % Solve least squares for coefficients
    [Q,R] = qr(Phi,0);
    coeffs = R \ (Q' * samples(:));

    % Reconstruct signal
    recon = real(Phi * coeffs);

    % Plot original vs reconstruction
    t = linspace(0, 2*pi, N);
    figure;
    plot(t, samples, 'LineWidth', 3, 'Color', 'g', 'DisplayName', 'Original');
    hold on;
    plot(t, recon, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Reconstruction');
    legend;
    title(sprintf('Reconstruction with single non-L2 basis: %s', metric_def_to_string(metric_def)));

    % Report error
    err = norm(samples(:) - recon(:)) / norm(samples(:));
    fprintf('Relative error with %s basis: %.6f\n', metric_def_to_string(metric_def), err);
end
