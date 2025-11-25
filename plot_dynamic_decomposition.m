function plot_dynamic_decomposition(coeffs, freqs, metrics, ax, n_samples)

    if ~exist('ax', 'var')
        figure;
        decomp_plot = subplot(1,1,1);
        hold(decomp_plot, 'on')
    else
        decomp_plot = ax;
        hold(decomp_plot, 'on')
    end
    
    if ~exist('n_samples', 'var')
        % oversample a bit so we can see the pretty shapes
        %n_samples = (round(max(freqs)/100) + 3) * 100;
        n_samples = max(freqs);
    end
        
    for i = 1:1:length(coeffs)
    
        single_signal = dynamic_basis_recomposition(coeffs(i), freqs(i), {metrics{i}}, n_samples);

        plot(decomp_plot, single_signal);


    end

    full_signal = dynamic_basis_recomposition(coeffs, freqs, metrics, n_samples);

    plot(decomp_plot, full_signal);

    hold(decomp_plot, 'off')
    
end