samples = randn(1,20);
%metrics = {make_weighted_p_metric(1,1,2), make_weighted_p_metric(10,1,2), make_weighted_p_metric(1,10,2)}
metrics = {make_weighted_p_metric(1,1,1), make_weighted_p_metric(1,1,2), make_weighted_p_metric(1,1,3), make_weighted_p_metric(1,1,4)}

figure;
samp_plot = subplot(2,2,1);
title(samp_plot, "Sample and Reconstructions")
hold(samp_plot, "on")
ratio_plot = subplot(2,2,2);
title(ratio_plot, "Ratio of Reconstructions over Samples")
hold(ratio_plot, "on")

% create vector of frequencies to plot over
freqs = ((0:length(samples)-1)*100)/length(samples);
mag_plot = subplot(2,2,3);
title(mag_plot, "P-DFT Magnitudes")
hold(mag_plot, "on")
phase_plot = subplot(2,2,4);
title(phase_plot, "P-DFT Phases")
hold(phase_plot, "on")

plot(samp_plot, samples, LineWidth=3, LineStyle="-", DisplayName="Original Samples")

for i = 1:1:(length(metrics))

    cur_metric = metrics{i};

    p_dft = alt_disc_fourier(samples, cur_metric);

    recon_samples_p = real(alt_inv_disc_fourier(p_dft, cur_metric));

    err = norm(samples - recon_samples_p);

    recon_legend_name = sprintf("Reconstructed with norm: %s , total error: %.3f", func2str(cur_metric), err);

    ratio = recon_samples_p./samples;

    plot(samp_plot, recon_samples_p, LineWidth=2, LineStyle="--", DisplayName=recon_legend_name)
    
    ratio_legend_name = sprintf("Ratio between reconstructed and original samples with norm: %s", func2str(cur_metric));

    plot(ratio_plot, ratio, LineWidth=2, LineStyle="--", DisplayName=ratio_legend_name);

    mag_legend_name = sprintf("Magnitude of DFT with norm: %s", func2str(cur_metric));

    plot(mag_plot, freqs, abs(p_dft), LineWidth=2, LineStyle="--", DisplayName=mag_legend_name);

    phase_legend_name = sprintf("Phase of DFT with norm: %s", func2str(cur_metric));

    % surpress teeny tiny coefficients
    surpressed_p_dft = p_dft;
    surpressed_p_dft(abs(p_dft)<1e-6) = 0;

    % and then extract the phases 
    phases = unwrap(angle(surpressed_p_dft));

    plot(phase_plot, freqs, phases, LineWidth=2, LineStyle="--", DisplayName=phase_legend_name);

end