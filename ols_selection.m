function [active, cA, infoA] = ols_selection(D, dict_info, x, percent_error_threshold, max_allowed_coeffs)
    if ~exist('max_allowed_coeffs','var'), max_allowed_coeffs = size(D,2); end
    x = x(:);
    r = x;
    active = [];
    cA = [];  % initialize coefficients
    threshold = percent_error_threshold * norm(x);

    while norm(r) > threshold && numel(active) < max_allowed_coeffs
        best_err = inf; best_idx = [];
        for k = 1:size(D,2)
            if ~isempty(active) && any(active == k), continue; end
            PhiA_test = D(:, [active, k]);
            [Q,R] = qr(PhiA_test, 0);
            c_test = R \ (Q' * x);
            r_test = x - PhiA_test * c_test;
            err = norm(r_test);
            if err < best_err
                best_err = err;
                best_idx = k;
            end
        end
        if isempty(best_idx), break; end
        active = [active, best_idx];
        PhiA = D(:, active);
        [Q,R] = qr(PhiA, 0);
        cA = R \ (Q' * x);
        r = x - PhiA * cA;
    end

    infoA = dict_info(active);
end

