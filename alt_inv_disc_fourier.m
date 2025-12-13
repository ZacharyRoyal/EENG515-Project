function recon = alt_inv_disc_fourier(coeffs, metric_def)
    Phi = build_alt_basis(N, metric_def);
    recon = Phi * coeffs;
end