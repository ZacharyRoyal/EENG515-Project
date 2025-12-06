metric_1_def = struct('W', [1 2; 0 2]);
metric_1_handle = make_weighted_p_metric_struct(metric_1_def);
metric_2_def = struct('W', [1 2; 0 2]);
metric_2_handle = make_weighted_p_metric_struct(metric_2_def);
plot_alt_unitball(metric_1_handle)
plot_alt_unitball(metric_2_handle)
frequency_1 = 2;
frequency_2 = 1;
wave_1_handle = @(t) alt_sin(frequency_1 * t, metric_1_handle)
wave_2_handle = @(t) alt_sin(frequency_2 * t, metric_2_handle)
integrand = @(theta) wave_1_handle(theta) .* wave_2_handle(theta)
inner_product = integral(integrand, 0, 2*pi)
