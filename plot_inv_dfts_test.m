% define signal here
step_size = 0.1;
time_series = 0:step_size:2*pi;
frequency = 4;

% random sample with varying smoothness
%random_smoothing = 2;
%samples = smoothed_random_samples(length(time_series), random_smoothing);

% triangle wave
%samples = tri_samples(time_series, frequency);

% sawtooth wave
%samples = saw_samples(time_series, frequency);

% square wave with intermediary flat zones at 0
%samples = square_samples(time_series, frequency);

% square wave which directly transitions from one plateu to the next
%samples = block_samples(time_series, frequency);

% funky shark-fin looking shapes
samples = shark_samples(time_series, frequency);

% sin wave of given metric
%samples = sin_samples(time_series, frequency, make_weighted_p_metric(1,1,0.5));
%samples = samples + sin_samples(time_series, frequency*3, make_weighted_p_metric(100,1,3));
%amples = samples + sin_samples(time_series, frequency*2, make_weighted_p_metric(1,1,2)) .* 0.2;

% add some random noise if you choose
%samples = add_random_noise(samples, 0.1, 5);

% increasingly square-ish using weighting of l2
%metric_defs = {struct('x_w', 1, 'y_w', 1), struct('x_w', 1/1.5, 'y_w', 1.5), struct('x_w', 1/2, 'y_w', 2), struct('x_w', 1/3, 'y_w', 3), struct('x_w', 1/4, 'y_w', 4)}

% increasingly impulse-like using weighting of l2
%metric_defs = {struct('x_w', 1, 'y_w', 1), struct('x_w', 1.5, 'y_w', 1/1.5), struct('x_w', 2, 'y_w', 1/2), struct('x_w', 3, 'y_w', 1/3), struct('x_w', 4, 'y_w', 1/4)}

% increasingly square-like sine waves using higher p-norm
metric_defs = {struct('p',0.5), struct('p',2), struct('p',5)};

% increasingly impulse-like sine waves using lower p-norm
%metric_defs = {struct('p',2), struct('p',1.5), struct('p',1), struct('p',0.5), struct('p',0.25)};

%full sweep from 0.25 to 10
%metric_defs = {struct('p',0.25), struct('p',0.5), struct('p',1), struct('p',2), struct('p',3), struct('p',4), struct('p',5), struct('p',6), struct('p',7), struct('p',8), struct('p',9), struct('p',10)};

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

for i = 1:1:(length(metric_defs))

    cur_metric = metric_defs{i};
    fprintf("%s\n", metric_def_to_string(metric_defs{i}));

    p_dft = alt_disc_fourier(samples, cur_metric);

    % 0 out some number of smallest coefficents
    %top_percent_to_keep = 0.50;
    %n = length(samples) - round(top_percent_to_keep*length(samples));
    %n = length(samples)-1; %for only most important componenet
    %[~, sorted_indexes] = sort(p_dft);
    %p_dft(sorted_indexes(1:n)) = 0;

    recon_samples_p = real(alt_inv_disc_fourier(p_dft, cur_metric));%real(alt_inv_disc_fourier(p_dft, cur_metric));

    err = norm(samples - recon_samples_p);

    recon_legend_name = sprintf("Reconstructed with p-norm: %.2f, total error: %.3f", cur_metric.p, err);

    ratio = samples./recon_samples_p;

    plot(samp_plot, recon_samples_p, LineWidth=2, LineStyle="--", DisplayName=recon_legend_name)
    
    ratio_legend_name = sprintf("Ratio between reconstructed and original samples with p-norm: %.2f", cur_metric.p);

    plot(ratio_plot, ratio, LineWidth=2, LineStyle="--", DisplayName=ratio_legend_name);

    mag_legend_name = sprintf("Magnitude of DFT with p-norm: %.2f", cur_metric.p);

    % surpress teeny tiny coefficients
    surpressed_p_dft = p_dft;
    surpressed_p_dft(abs(p_dft)<1e-6) = 0;

    plot(mag_plot, freqs, abs(surpressed_p_dft), LineWidth=2, LineStyle="--", DisplayName=mag_legend_name);

    phase_legend_name = sprintf("Phase of DFT with p-norm: %.2f", metric_defs{i}.p);

    % and then extract the phases 
    phases = unwrap(angle(surpressed_p_dft));

    plot(phase_plot, freqs, phases, LineWidth=2, LineStyle="--", DisplayName=phase_legend_name);

end