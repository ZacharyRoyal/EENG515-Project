function [coeffs, freqs, norms] = enhanced_dynamic_basis_decomposition(samples, percent_error_threshold)

    % initialize return variables, these will be grown into arrays later
    coeffs = [];
    freqs = [];
    norms = [];

    % define all the metrics we will check over
    %p_sweep = 0.2:0.2:5
    %p_sweep = 0.1:0.1:5;
    %p_sweep = [0.8, 1:0.5:6];
    %p_sweep = 0.5:0.5:6; % current champion for convergence robustness and term count
    %p_sweep = [0.05, 0.1, 0.2, 0.4, 0.8, 1, 1.5, 2, 3, 4, 5];
    %p_sweep = [0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5 4, 4.5, 5, 10];
    p_sweep = [2]; %for comparing to base DFT
    %p_sweep = [10]
    %p_sweep = 0.5:0.25:5;
    %p_sweep = [1,2,3];

    metric_def_sweep = cell(size(p_sweep));

    for i = 1:1:length(p_sweep)
        metric_def_sweep{i} = struct('p', p_sweep(i));
    end

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

    while (norm(error_signal) > energy_threshold)
        
        % now we loop over each metric, see which one has the best 1-term
        % approximation, and then subtract that approximation from the
        % signal, and rinse repeat until the total energy drops below a
        % threshold

        best_error = 10^10;
        best_p_index = 1;
        best_coeff = 0;
        best_freq = 0;

        for i = 1:1:length(metric_def_sweep)
            
            current_metric = make_weighted_p_metric_struct(metric_def_sweep{i});

            % decompose the signal
            p_dft = alt_disc_fourier(samples, current_metric);

            for j = 1:1:length(p_dft)
                current_term = p_dft(j);
                current_freq = j;

                candidate_coeffs = [coeffs, current_term];
                candidate_freqs = [freqs, current_freq];
                candidate_norms = [norms, p_sweep(i)];

                candidate_signal = dynamic_basis_recomposition(candidate_coeffs, candidate_freqs, candidate_norms, length(samples));
                candidate_error = norm(samples - candidate_signal);

                if candidate_error < best_error
                    best_error = candidate_error;
                    best_p_index = i;
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
        norms(term_index) = p_sweep(best_p_index);

        current_signal = dynamic_basis_recomposition(coeffs, freqs, norms, length(samples));

        error_signal = samples - current_signal;

        error_norms(end+1) = norm(error_signal);
        threshold_line = ones(size(error_norms)).*(energy_threshold);

        plot(error_plot, error_norms, LineWidth=1.5, LineStyle='--');
        hold(error_plot, 'on')
        plot(error_plot, threshold_line, LineWidth=1)
        hold(error_plot, 'off')

        pause(0.01)

        term_index = term_index + 1;

    end

end