% Instead of solving the normal equations directly by constructing the Gram
% matrix G = \Phi^{T}*\Phi
% Use QR factorization on \Phi, so \Phi = QR, G is never constructed
% explicitly. This increases numerical stability
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
        PhiA = D(:, active);    % Subset of dictionary atoms (Phi)
        [Q,R] = qr(PhiA, 0);    % Factor Phi = Q R
        cA = R \ (Q' * x);      % Solve R*c=Q'*x -> c=R^{-1}*Q'*x
        r = x - PhiA * cA;      % Residual after reconstruction

        % If you instead want to solve normal equations directly
        % PhiA = D(:, active);        % Subset of dictionary atoms, active (Phi)
        % G = PhiA' * PhiA;           % Gram matrix
        % b = PhiA' * x;              % right-hand side
        % cA = b / G;                 % solve b = G^(-1)*c
        % r = x - PhiA * cA;          % residual

    end
end
