function plot_dynamic_decomposition(coeffs, freqs, metrics, ax)

    if ~exist('ax', 'var')
        figure;
        decomp_plot = subplot(1,1,1);
        hold(decomp_plot, 'on')
    else
        decomp_plot = ax;
        hold(decomp_plot, 'on')
    end
    
    % oversample a bit so we can see the pretty shapes
    n_samples = (round(max(freqs)/100) + 3) * 100;
    
    for i = 1:1:length(coeffs)
    
        single_signal = dynamic_basis_recomposition(coeffs(i), freqs(i), {metrics{i}}, n_samples);

        plot(decomp_plot, single_signal);


    end

    hold(decomp_plot, 'off')
    
end