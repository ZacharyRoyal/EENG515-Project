function metric_handle = make_weighted_p_metric_struct(metric_def)

    % default to unweighted l2
    x_w = 1;
    y_w = 1;
    p = 2;
    
    if isfield(metric_def, 'x_w')
        x_w = metric_def.x_w;
    end

    if isfield(metric_def, 'y_w')
        y_w = metric_def.y_w;
    end

    if isfield(metric_def, 'p')
        p = metric_def.p;
    end

    metric_handle = @(x, y) ( ...
        (abs(x)*x_w).^p + (abs(y)*y_w).^p...
        ).^(1/p);
end