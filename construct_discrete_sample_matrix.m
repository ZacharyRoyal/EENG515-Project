function Phi = construct_discrete_sample_matrix(N, metric_def)

    Phi = zeros(N);

    metric_handle = make_weighted_p_metric_struct(metric_def);

    for k = 0:1:N-1 % frequency index

        for n = 0:1:N-1 % sample index counter 

            power_term = 2 * pi * (k/N) * n;

            current_inc = alt_e_power_i_x(power_term, metric_handle);
            Phi(k+1, n+1) = Phi(k+1, n+1) + current_inc;

        end

    end

end