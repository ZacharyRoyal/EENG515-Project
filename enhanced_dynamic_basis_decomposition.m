function [coeffs, freqs, metrics] = enhanced_dynamic_basis_decomposition(samples, metric_def_list, percent_error_threshold, max_allowed_coeffs)

    if ~exist('max_allowed_coeffs','var'), max_allowed_coeffs = numel(samples); end
    N = numel(samples);

    % Build a dictionary of atoms (columns) across metrics and frequencies
    all_freqs = 0:(N-1);
    dict_cols = {};
    dict_info = {};
    for m = 1:numel(metric_def_list)
        Phi_m = build_alt_basis(N, all_freqs, metric_def_list{m});
        % normalize for fair correlation
        Phi_m = bsxfun(@rdivide, Phi_m, sqrt(sum(abs(Phi_m).^2,1)));
        for j = 1:numel(all_freqs)
            dict_cols{end+1} = Phi_m(:,j);
            dict_info{end+1} = struct('freq', all_freqs(j), 'metric', metric_def_list{m});
        end
    end
    D = cell2mat(dict_cols);  % N x K

    r = samples(:);
    active = [];
    coeffs = [];
    freqs = [];
    metrics = {};
    threshold = percent_error_threshold * norm(samples);

    while norm(r) > threshold && numel(active) < max_allowed_coeffs
        % Select atom with max absolute correlation
        corr = abs(D' * r);
        [~, idx] = max(corr);
        active = unique([active, idx]);  % avoid duplicates

        PhiA = D(:, active);
        % Refit coefficients via QR
        [Q,R] = qr(PhiA,0);
        cA = R \ (Q' * samples(:));
        xhat = PhiA * cA;
        r = samples(:) - xhat;

        % Update outputs
        coeffs = cA;
        freqs = cellfun(@(s) s.freq, dict_info(active));
        metrics = cellfun(@(s) s.metric, dict_info(active), 'UniformOutput', false);
    end
end
