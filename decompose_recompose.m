% define signal here
step_size = 0.1;
time_series = 0:step_size:2*pi;
frequency = 2;

% random sample with varying smoothness
%random_smoothing = 1;
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
%samples = sin_samples(time_series, frequency, make_weighted_p_metric(1,1,2));

% add some random noise if you choose
%samples = add_random_noise(samples, 0.1, 0);

percent_error_threshold = 0.1;

figure;
recomp_compare_plot = subplot(1,1,1);
plot(recomp_compare_plot, samples, LineWidth=4, Color='g')
hold(recomp_compare_plot, 'on')

figure;
decomp_compare_plot_1 = subplot(3,1,1);
decomp_compare_plot_2 = subplot(3,1,2);
decomp_compare_plot_3 = subplot(3,1,3);

% first we evaluate the standard l2 decomposition
l2_coeffs_base = fft(samples);
percent_error = 10e10;
samples_to_zero = length(l2_coeffs_base);
% just creating this variable out of the loop so it can be used below it
% later, it's value will be overwritten in the loop
recon_samples_l2 = real(ifft(l2_coeffs_base));
while (percent_error > percent_error_threshold)

    % increment our allowed samples
    samples_to_zero = max(0, samples_to_zero - 1);

    % sort coeffs by magnitude
    [~, sorted_indexes] = sort(abs(l2_coeffs_base));

    l2_coeffs_cut = l2_coeffs_base;

    % 0 out some number of the smallest ones
    l2_coeffs_cut((sorted_indexes(1:samples_to_zero))) = 0;

    % reconstruct the signal with the remaining ones
    recon_samples_l2 = real(ifft(l2_coeffs_cut));

    % and see how close we are
    percent_error = norm(samples - recon_samples_l2)/norm(samples);

end

l2_terms_to_approx = length(l2_coeffs_base) - samples_to_zero
percent_error = norm(samples - recon_samples_l2)/norm(samples)
plot(recomp_compare_plot, recon_samples_l2, LineWidth=2, LineStyle="--")
l2_list = {};
for i = 1:1:length(l2_coeffs_cut)
    l2_list{i} = struct('p', 2);
end
plot_dynamic_decomposition(l2_coeffs_cut, 1:1:length(l2_coeffs_cut), l2_list, decomp_compare_plot_1, length(samples));

% define all the metrics we will check over
%p_sweep = 0.2:0.2:5
%p_sweep = 0.1:0.1:5;
%p_sweep = [0.8, 1:0.5:6];
%p_sweep = 0.5:0.5:6; % current champion for convergence robustness and term count
%p_sweep = [0.05, 0.1, 0.2, 0.4, 0.8, 1, 1.5, 2, 3, 4, 5];
%p_sweep = [0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5 4, 4.5, 5, 10];
%p_sweep = [2]; %for comparing to base DFT
%p_sweep = [10];
%p_sweep = 0.5:0.25:5;
%p_sweep = [1,2,3];
%p_sweep = 1:0.5:6;

%metric_def_sweep = cell(size(p_sweep));

%for i = 1:1:length(p_sweep)
%    metric_def_sweep{i} = struct('p', p_sweep(i));
%end

%w_sweep = 0.5:0.5:5; %current champion 
%w_sweep = 0.1:0.1:10; % danger zone
%w_sweep = [0.1, 0.5, 1, 2, 4, 8];
w_sweep = [0.01, 1, 100]; % new champion?
%w_sweep = [100];
%w_sweep = [1]; % also for comparing to base DFT
metric_def_sweep = cell(size(w_sweep));
for i = 1:1:length(w_sweep)
    metric_def_sweep{i} = struct('x_w', w_sweep(i));
end

% then the first iteraiton of the dynamic basis decomposition
[coeffs, freqs, metrics] = dynamic_basis_decomposition(samples, metric_def_sweep, percent_error_threshold);
dyn_terms_to_approx = length(coeffs)
recon_samples_dyn = dynamic_basis_recomposition(coeffs, freqs, metrics, length(samples));
percent_error = norm(samples - recon_samples_dyn)/norm(samples)

plot(recomp_compare_plot, recon_samples_dyn, LineWidth=2, LineStyle="--")
plot_dynamic_decomposition(coeffs, freqs, metrics, decomp_compare_plot_2, length(samples));

% and finally the enhanced dynamic basis decomposition
[coeffs, freqs, metrics] = enhanced_dynamic_basis_decomposition(samples, metric_def_sweep, percent_error_threshold);
enhanced_dyn_terms_to_approx = length(coeffs)
recon_samples_enhanced_dyn = dynamic_basis_recomposition(coeffs, freqs, metrics, length(samples));
percent_error = norm(samples - recon_samples_enhanced_dyn)/norm(samples)

plot(recomp_compare_plot, recon_samples_enhanced_dyn, LineWidth=2, LineStyle="--")
plot_dynamic_decomposition(coeffs, freqs, metrics, decomp_compare_plot_3, length(samples));