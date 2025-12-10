function Phi = build_alt_basis(N, freqs, metric_def)
    metric_handle = make_weighted_p_metric_struct(metric_def);
    Phi = zeros(N, numel(freqs));
    for j = 1:numel(freqs)
        f = freqs(j);
        for n = 0:N-1
            theta = 2*pi*(f/N)*n;
            Phi(n+1, j) = alt_e_power_i_x(theta, metric_handle);
        end
    end
end
