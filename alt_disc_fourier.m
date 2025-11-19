function coeffs = alt_disc_fourier(samples, metric_handle)

    N = size(samples, 2);
    coeffs = zeros(1, N);

    for k = 0:1:N-1
        
        for n = 0:1:N-1

            power_term = -2 * pi * (k/N) * n;

            alpha = alt_alpha(power_term, metric_handle);
            alpha_correction = 1/(alpha^2);
                
            sample_inc = samples(n+1) * alt_e_power_i_x(power_term, metric_handle);
            sample_inc = sample_inc * alpha_correction;

            coeffs(k+1) = coeffs(k+1) + sample_inc;
        
        end

    end

end