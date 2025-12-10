function test_signals()

    % ----- Parameters -----
    N = 100;                 % number of samples
    t = linspace(0, 2*pi, N);
    percent_error_threshold = 0.1;

    % ----- Define signals -----
    signals = struct();
    signals.sine      = sin(4*t);
    signals.square    = square(4*t);
    signals.sawtooth  = sawtooth(4*t);
    signals.triangle  = sawtooth(4*t, 0.5);
    signals.noise     = randn(1,N);
    signals.shark     = shark_samples(t, 4); % your custom function

    % ----- Metrics sweep -----
    metric_def_list = {
        struct('p', 1), struct('p', 2), struct('p', 3), struct('W', [0.01 0; 0 1]), struct('W', [100 0; 0 1]), struct('W', [1 0.5; 0 1]), struct('W', [1 -0.5; 0 1])
    };

    % ----- Loop over signals -----
    signal_names = fieldnames(signals);
    summary = {}; % collect results

    for s = 1:numel(signal_names)
        name = signal_names{s};
        x = signals.(name)(:);

        % Build dictionary
        freqs = 0:(N-1);
        [D, dict_info] = build_alt_dictionary(N, freqs, metric_def_list);
        D = D ./ vecnorm(D); % normalize columns

        % --- Run OMP ---
        [active_idx_omp, coeffs_omp] = omp_qr(D, x, percent_error_threshold, N);
        PhiA_omp = D(:, active_idx_omp);
        recon_omp = real(PhiA_omp * coeffs_omp);
        omp_terms = numel(active_idx_omp);
        omp_error = norm(x - recon_omp) / norm(x);
        recon_omp = recon_omp(:);

           
        % --- Run OLS ---
        [active_idx_ols, coeffs_ols, ~] = ols_selection(D, dict_info, x, percent_error_threshold, N);
        PhiA_ols = D(:, active_idx_ols);
        recon_ols = real(PhiA_ols * coeffs_ols);
        ols_terms = numel(active_idx_ols);
        ols_error = norm(x - recon_ols) / norm(x);
        recon_ols = recon_ols(:);

        % --- FFT baseline ---
        fft_coeffs = fft(x);
        percent_error = inf;
        samples_to_zero = length(fft_coeffs);
        recon_fft = real(ifft(fft_coeffs));
        recon_fft = recon_fft(:);
        while percent_error > percent_error_threshold
            samples_to_zero = max(0, samples_to_zero - 1);
            [~, sorted_idx] = sort(abs(fft_coeffs));
            fft_cut = fft_coeffs;
            fft_cut(sorted_idx(1:samples_to_zero)) = 0;
            recon_fft = real(ifft(fft_cut));
            percent_error = norm(x - recon_fft) / norm(x);
        end
        l2_terms = length(fft_coeffs) - samples_to_zero;
        l2_error = percent_error;

        % --- Plot ---
        figure('Name', ['Signal: ', name]);
        plot(t, x, 'LineWidth', 3, 'Color', 'g', 'DisplayName', 'Original');
        hold on;
        plot(t, recon_fft, 'LineWidth', 2, 'LineStyle', '--', ...
            'DisplayName', sprintf('FFT (%d terms, err=%.3f)', l2_terms, l2_error));
        plot(t, recon_omp, 'LineWidth', 2, 'LineStyle', '--', ...
            'DisplayName', sprintf('OMP (%d terms, err=%.3f)', omp_terms, omp_error));
        plot(t, recon_ols, 'LineWidth', 2, 'LineStyle', '--', ...
     'DisplayName', sprintf('OLS (%d terms, err=%.3f)', ols_terms, ols_error));
        legend;
        title(['OMP vs OLS vs FFT reconstruction for ', name]);

        % --- Save summary ---
        summary{end+1} = struct('signal', name, ...
                                'fft_terms', l2_terms, 'fft_error', l2_error, ...
                                'omp_terms', omp_terms, 'omp_error', omp_error, ...
                                'ols_terms', ols_terms, 'ols_error', ols_error);
    end

    % ----- Print summary table -----
    fprintf("\nSummary across signals:\n");
    fprintf("%-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s\n", ...
        'Signal', 'FFT terms', 'FFT error', 'OMP terms', 'OMP error', 'OLS terms', 'OLS error');
    fprintf("-------------------------------------------------------------------------\n");
    for i = 1:numel(summary)
        s = summary{i};
        fprintf("%-10s | %-10d | %-10.4f | %-10d | %-10.4f | %-10d | %-10.4f\n", ...
            s.signal, s.fft_terms, s.fft_error, s.omp_terms, s.omp_error, s.ols_terms, s.ols_error);
    end
end

