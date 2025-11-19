% define signal here
time_series = 0:0.1:2*pi;
frequency = 4;

% random sample with varying smoothness
random_smoothing = 15;
samples = smoothed_random_samples(length(time_series), random_smoothing);

% triangle wave
%samples = tri_samples(time_series, frequency);

% sawtooth wave
%samples = saw_samples(time_series, frequency);

% square wave
%samples = square_samples(time_series, frequency);

% sin wave of given metric
%samples = sin_samples(time_series, frequency, make_weighted_p_metric(5,1,3));

percent_error_threshold = 0.01;

figure;
recomp_compare_plot = subplot(1,1,1);
%norm_histogram_plot = subplot(2,1,2);
plot(recomp_compare_plot, samples, LineWidth=4, Color='g')
hold(recomp_compare_plot, 'on')
%hold(norm_histogram_plot, 'on')

% first we evaluate the standard l2 decomposition
[coeffs, freqs, norms] = l2_decomposition(samples, percent_error_threshold);
l2_terms_to_approx = length(coeffs)
recon_samples_l2 = dynamic_basis_recomposition(coeffs, freqs, norms, length(samples));
percent_error = norm(samples - recon_samples_l2)/norm(samples)

plot(recomp_compare_plot, recon_samples_l2, LineWidth=2, LineStyle="--")

% then the first iteraiton of the dynamic basis decomposition
%[coeffs, freqs, norms] = dynamic_basis_decomposition(samples, percent_error_threshold);
%dyn_terms_to_approx = length(coeffs)
%recon_samples_dyn = dynamic_basis_recomposition(coeffs, freqs, norms, length(samples));
%percent_error = norm(samples - recon_samples_dyn)/norm(samples)

%plot(recomp_compare_plot, recon_samples_dyn, LineWidth=2, LineStyle="--")
%plot(norm_histogram_plot, sort(norms));

% and finally the enhanced dynamic basis decomposition
[coeffs, freqs, norms] = enhanced_dynamic_basis_decomposition(samples, percent_error_threshold);
enhanced_dyn_terms_to_approx = length(coeffs)
recon_samples_enhanced_dyn = dynamic_basis_recomposition(coeffs, freqs, norms, length(samples));
percent_error = norm(samples - recon_samples_enhanced_dyn)/norm(samples)

plot(recomp_compare_plot, recon_samples_enhanced_dyn, LineWidth=2, LineStyle="--")

%plot(norm_histogram_plot, sort(norms));