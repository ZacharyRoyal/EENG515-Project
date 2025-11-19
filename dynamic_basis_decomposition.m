function [coeffs, freqs, norms] = dynamic_basis_decomposition(samples)

    % initialize return variables, these will be grown into arrays later
    coeffs = 0;
    freqs = 0;
    norms = 0;

    % first define all the metrics we will check over
    p_sweep = 0.25:0.25:5;
    %p_sweep = 1:1:5;
    metric_def_sweep = cell(size(p_sweep));

    for i = 1:1:length(p_sweep)
        metric_def_sweep{i} = struct('p', p_sweep(i));
    end

    % break into local variable for convenience
    signal = samples;
    energy_threshold = 0.05*norm(signal); %threshold for good enough
    % currently 5% of original signal energy

    term_index = 1;

    % debug variable to allow plotting of signal norm over time 
    signal_norms = [norm(signal)];

    while (norm(signal) > energy_threshold)
        
        % now we loop over each metric, see which one has the best 1-term
        % approximation, and then subtract that approximation from the
        % signal, and rinse repeat until the total energy drops below a
        % threshold

        best_error = 10^10
        best_p_index = 1;
        best_coeff = 0;
        best_freq = 0;

        for i = 1:1:length(metric_def_sweep)
            
            current_metric = make_weighted_p_metric_struct(metric_def_sweep{i});

            % decompose the signal
            p_dft = alt_disc_fourier(signal, current_metric)

            % sort terms by magnitude
            [sorted_terms, sorted_indexes] = sort(abs(p_dft))

            max_index = sorted_indexes(end)

            dominant_term = p_dft(max_index)

            % zero-out non-dominant terms
            p_dft(sorted_indexes(1:end-1)) = 0

            reconstructed_signal = real(alt_inv_disc_fourier(p_dft, current_metric))

            error = norm(signal - reconstructed_signal)

            if error < best_error
                best_error = error
                best_p_index = i;
                best_coeff = dominant_term;
                best_freq = max_index;

                best_p = metric_def_sweep{i}.p
                plot(signal)
                hold on
                plot(reconstructed_signal)
                hold off
                pause(0.1)

            end

        end
        
        coeffs(term_index) = best_coeff
        freqs(term_index) = best_freq
        norms(term_index) = p_sweep(best_p_index)

        best_dft = zeros(1,length(signal))
        best_dft(best_freq) = best_coeff

        best_reconstructed_signal = real(alt_inv_disc_fourier(best_dft, make_weighted_p_metric_struct(metric_def_sweep{best_p_index})));

        plot(signal)
        hold on
        plot(best_reconstructed_signal)

        new_signal = signal - best_reconstructed_signal;
        plot(new_signal)
        hold off

        signal_norms(end+1) = norm(new_signal)

        pause(0.1)

        term_index = term_index + 1;
        signal = new_signal

    end

end