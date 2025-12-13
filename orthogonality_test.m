metric_def = struct('W', [1 2; 0 2]);
metric_handle = make_weighted_p_metric_struct(metric_def);
plot_alt_unitball(metric_handle)
frequency_1 = 5;
frequency_2 = 6;
wave_1_handle = @(t) alt_sin(frequency_1 * t, metric_handle)
wave_2_handle = @(t) alt_sin(frequency_2 * t, metric_handle)
integrand = @(theta) wave_1_handle(theta) .* wave_2_handle(theta)
inner_product = integral(integrand, 0, 2*pi)
