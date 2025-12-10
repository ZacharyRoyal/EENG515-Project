function coeffs = alt_disc_fourier(samples, metric_def)
    N = numel(samples);
    Phi = build_alt_basis(N, freqs, metric_def);
    coeffs = analyze_ls(Phi, samples);
end