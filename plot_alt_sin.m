function plot_alt_sin(theta_vals, metric_handle)

    theta_count = max(size(theta_vals));

    vals = zeros(theta_count, 1);

    for i = 1:1:theta_count
        
        vals(i) = alt_sin(theta_vals(i), metric_handle);

    end

    figure;

    plot(theta_vals, vals);

end