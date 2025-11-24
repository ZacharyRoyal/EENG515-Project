function [coeffs, freqs, metrics] = enhanced_dynamic_basis_decomposition(samples, metric_def_list, percent_error_threshold, max_allowed_coeffs)

    % initialize return variables, these will be grown into arrays later
    coeffs = [];
    freqs = [];
    metrics = {};

    energy_threshold = percent_error_threshold*norm(samples); %threshold for good enough
    % basically X% of original signal energy

    term_index = 1;

    current_signal = zeros(size(samples));
    error_signal = samples - current_signal;

    % debug variable to allow plotting of error norm over time 
    error_norms = [norm(error_signal)];

    figure;
    signal_plot = subplot(2,1,1);
    error_plot =subplot(2,1,2);
    
    if ~exist('max_allowed_coeffs', 'var')
        max_allowed_coeffs = length(samples);
    end

    coeff_count = 0;

    while (norm(error_signal) > energy_threshold && coeff_count <= max_allowed_coeffs)
        
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
            p_dft = alt_disc_fourier(samples, current_metric);

            for j = 1:1:length(p_dft)
                current_term = p_dft(j);
                current_freq = j;

                candidate_coeffs = [coeffs, current_term];
                candidate_freqs = [freqs, current_freq];
                candidate_metrics = metrics;
                candidate_metrics{end+1} = current_metric;

                candidate_signal = dynamic_basis_recomposition(candidate_coeffs, candidate_freqs, candidate_metrics, length(samples));
                candidate_error = norm(samples - candidate_signal);

                if candidate_error < best_error
                    best_error = candidate_error;
                    best_metric = current_metric;
                    best_coeff = current_term;
                    best_freq = current_freq;
    
                    plot(signal_plot, samples, LineWidth=3)
                    hold(signal_plot, 'on')
                    plot(signal_plot, candidate_signal, LineWidth=2, LineStyle="--")
                    hold(signal_plot, 'off')
                    pause(0.01)
    
                end

            end

        end
        
        coeffs(term_index) = best_coeff;
        freqs(term_index) = best_freq;
        metrics{term_index} = best_metric;

        current_signal = dynamic_basis_recomposition(coeffs, freqs, metrics, length(samples));

        error_signal = samples - current_signal;

        error_norms(end+1) = norm(error_signal);
        threshold_line = ones(size(error_norms)).*(energy_threshold);

        plot(error_plot, error_norms, LineWidth=1.5, LineStyle='--');
        hold(error_plot, 'on')
        plot(error_plot, threshold_line, LineWidth=1)
        hold(error_plot, 'off')

        pause(0.01)

        term_index = term_index + 1;
        coeff_count = coeff_count + 1;

    end

end