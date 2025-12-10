function [active, cA] = omp_qr(D, x, percent_error_threshold, max_allowed_coeffs)
    if ~exist('max_allowed_coeffs', 'var')
        max_allowed_coeffs = size(D,2); 
    end
    x = x(:);
    r = x;
    active = [];
    threshold = percent_error_threshold * norm(x);

    while norm(r) > threshold && numel(active) < max_allowed_coeffs
        % Select atom with max absolute correlation
        corr = abs(D' * r);
        [~, idx] = max(corr);
        if ~isempty(active) && any(active == idx)
            % If duplicated, break to avoid stalling
            break;
        end
        active = [active, idx];
        % Refit coefficients on active set using QR
        PhiA = D(:, active);
        [Q,R] = qr(PhiA, 0);
        cA = R \ (Q' * x);
        r = x - PhiA * cA;
    end
end
