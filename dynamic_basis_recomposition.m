function samples = dynamic_basis_recomposition(coeffs, freqs, norms, n_samples)

    % identify highest frequnecy
    max_freq = max(freqs);

    if exist('n_samples', 'var') == 1
        max_freq = n_samples;
    end

    samples = zeros(1,max_freq);

    for i = 1:1:length(coeffs)
        current_coeff = coeffs(i);
        current_freq = freqs(i);
        current_norm = norms(i);

        current_dft = zeros(1,max_freq);
        current_dft(current_freq) = current_coeff;

        currrent_signal = real(alt_inv_disc_fourier(current_dft, make_weighted_p_metric_struct(struct('p', current_norm))));

        %plot(samples)
        %hold on
        %plot(currrent_signal)        

        samples = samples + currrent_signal;
        %plot(samples)
        %hold off
        %pause(0.01)

    end

end