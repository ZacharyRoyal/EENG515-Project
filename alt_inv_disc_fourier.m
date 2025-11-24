function r_samps = alt_inv_disc_fourier(coeffs, metric_def)

    N = size(coeffs, 2);
    r_samps = zeros(1, N);

    metric_handle = make_weighted_p_metric_struct(metric_def);

    for n = 0:1:N-1
        
        for k = 0:1:N-1

            if coeffs(k+1) == 0
                continue
            end

            power_term = 2 * pi * (k/N) * n;

            coeff_inc = coeffs(k+1) * alt_e_power_i_x(power_term, metric_handle);

            r_samps(n+1) = r_samps(n+1) + coeff_inc;
        
        end

    end

    r_samps = r_samps/N;

end