function omp_decomp_recomp()

    % ----- Signal -----
    step_size = 0.1;
    t = 0:step_size:2*pi;
    frequency = 4;

    % Shark + flipped shark
    samples = shark_samples(t, frequency) + 1;
    % samples = samples + flip(samples);
    samples = samples(:); % column vector

    percent_error_threshold = 0.1;

    % ----- Metrics sweep (same spirit as your script) -----
    metric_def_list = {
        struct('p', 1), struct('p', 2), struct('p', 3), ...
        struct('W', [0.01 0; 0 1], 'p', 2), ...
        struct('W', [100 0; 0 1], 'p', 2), ...
        struct('W', [1 0.5; 0 1], 'p', 2), ...
        struct('W', [1 -0.5; 0 1], 'p', 2)
    };

    % ----- FFT baseline with coefficient pruning -----
    l2_coeffs = fft(samples);
    percent_error = inf;
    samples_to_zero = length(l2_coeffs);
    recon_l2 = real(ifft(l2_coeffs));
    while percent_error > percent_error_threshold
        samples_to_zero = max(0, samples_to_zero - 1);
        [~, sorted_indexes] = sort(abs(l2_coeffs));
        l2_cut = l2_coeffs;
        l2_cut(sorted_indexes(1:samples_to_zero)) = 0;
        recon_l2 = real(ifft(l2_cut));
        percent_error = norm(samples - recon_l2) / norm(samples);
    end
    l2_terms_to_approx = length(l2_coeffs) - samples_to_zero;

    % ----- Build dictionary across metrics and frequencies -----
    N = numel(samples);
    [D, dict_info] = build_alt_dictionary(N, metric_def_list); % N x K

    % ----- Run OMP with QR refit -----
    max_allowed_coeffs = N; % or smaller cap
    [active_idx, coeffs_omp] = omp_qr(D, samples, percent_error_threshold, max_allowed_coeffs);

    PhiA = D(:, active_idx);
    recon_omp = real(PhiA * coeffs_omp);
    dyn_terms_to_approx = numel(active_idx);
    percent_error_omp = norm(samples - recon_omp) / norm(samples);

    % ----- Report components used -----
    fprintf("Components used (OMP):\n");
    for k = 1:numel(active_idx)
        info = dict_info{active_idx(k)};
        fprintf("\tcoeff=%.3f%+ .3fi, freq=%d, %s\n", real(coeffs_omp(k)), imag(coeffs_omp(k)), info.freq, metric_def_to_string(info.metric));
    end

    % ----- Plots -----
    figure;
    recomp_compare_plot = subplot(1,1,1);
    plot(recomp_compare_plot, samples, 'LineWidth', 4, 'Color', 'g', 'DisplayName', 'Original');
    hold(recomp_compare_plot, 'on');
    plot(recomp_compare_plot, recon_l2, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', sprintf('FFT (%d terms, err=%.3f)', l2_terms_to_approx, norm(samples - recon_l2) / norm(samples)));
    plot(recomp_compare_plot, recon_omp, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', sprintf('OMP (%d terms, err=%.3f)', dyn_terms_to_approx, percent_error_omp));
    legend(recomp_compare_plot);
    title(recomp_compare_plot, 'Reconstruction comparison');

    % ----- Print summary -----
    fprintf("\nSummary:\n");
    fprintf("FFT terms: %d, error: %.4f\n", l2_terms_to_approx, norm(samples - recon_l2) / norm(samples));
    fprintf("OMP  terms: %d, error: %.4f\n", dyn_terms_to_approx, percent_error_omp);

end
