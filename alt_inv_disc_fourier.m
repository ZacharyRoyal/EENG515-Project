function r_samps = alt_inv_disc_fourier(coeffs, metric_handle)

    N = size(coeffs, 2);
    r_samps = zeros(1, N);

    for n = 0:1:N-1
        
        for k = 0:1:N-1

            power_term = 2 * pi * (k/N) * n;

            alpha = alt_alpha(power_term, metric_handle);
            alpha_correction = -(alpha)^2;

            coeff_inc = coeffs(k+1) * alt_e_power_i_x(power_term, metric_handle);
            %coeff_inc = coeffs(k+1) * alt_e_power_i_x(power_term * alpha_correction, metric_handle);

            r_samps(n+1) = r_samps(n+1) + coeff_inc;
        
        end

    end

    r_samps = r_samps/N;

end