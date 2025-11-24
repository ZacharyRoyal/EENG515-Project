function [coeffs, freqs, metrics] = dynamic_basis_decomposition(samples, metric_def_list, percent_error_threshold, max_allowed_coeffs)

    % initialize return variables, these will be grown into arrays later
    coeffs = [];
    freqs = [];
    metrics = {};

    % break into local variable for convenience
    signal = samples;
    energy_threshold = percent_error_threshold*norm(signal); %threshold for good enough
    % basically X% of original signal energy

    term_index = 1;

    % debug variable to allow plotting of signal norm over time 
    signal_norms = [norm(signal)];

    if ~exist('max_allowed_coeffs', 'var')
        max_allowed_coeffs = length(samples);
    end

    coeff_count = 0;

    figure;

    while (norm(signal) > energy_threshold && coeff_count <= max_allowed_coeffs)
        
        % now we loop over each metric, see which one has the best 1-term
        % approximation, and then subtract that approximation from the
        % signal, and rinse repeat until the total energy drops below a
        % threshold

        best_error = 10^10;
        best_metric = metric_def_list{1};
        best_coeff = 0;
        best_freq = 0;

        for i = 1:1:length(metric_def_list)
            
            current_metric = metric_def_list{i};

            % decompose the signal
            p_dft = alt_disc_fourier(signal, current_metric);

            % sort terms by magnitude
            [sorted_terms, sorted_indexes] = sort(abs(p_dft));

            max_index = sorted_indexes(end);

            dominant_term = p_dft(max_index);

            % zero-out non-dominant terms
            p_dft(sorted_indexes(1:end-1)) = 0;
            
            % old way just trying for best 1-term error
            reconstructed_signal = real(alt_inv_disc_fourier(p_dft, current_metric));

            error = norm(signal - reconstructed_signal);

            if error < best_error
                best_error = error;
                best_metric = current_metric;
                best_coeff = dominant_term;
                best_freq = max_index;

                plot(signal, LineWidth=3)
                hold on
                plot(reconstructed_signal, LineWidth=2, LineStyle="--")
                hold off
                pause(0.01)

            end

        end
        
        coeffs(term_index) = best_coeff;
        freqs(term_index) = best_freq;
        metrics{term_index} = best_metric;

        best_dft = zeros(1,length(signal));
        best_dft(best_freq) = best_coeff;

        best_reconstructed_signal = real(alt_inv_disc_fourier(best_dft, best_metric));

        new_signal = signal - best_reconstructed_signal;

        signal_norms(end+1) = norm(new_signal);

        pause(0.01)

        term_index = term_index + 1;
        signal = new_signal;
        coeff_count = coeff_count + 1;

    end

end