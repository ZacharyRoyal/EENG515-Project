function str = metric_def_to_string(metric_def)

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

    W = [x_w 0; 0 y_w];
    if isfield(metric_def, 'W')
        W = metric_def.W;
    end

    str = sprintf("p: %.2f, W: [%.2f, %.2f; %.2f, %.2f]", p, W(1,1), W(1,2), W(2,1), W(2,2));

end