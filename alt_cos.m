function c = alt_cos(theta, metric_handle)
    c = cos(theta) .* alt_alpha(theta, metric_handle);
end