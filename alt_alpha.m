function alpha = alt_alpha(theta, metric_handle)
    
    base_x = cos(theta);
    base_y = sin(theta);

    alpha = ones(size(theta))./metric_handle(base_x, base_y);

end