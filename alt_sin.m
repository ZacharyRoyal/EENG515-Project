function s = alt_sin(theta, metric_handle)
  
    s = (sin(theta).*ones(size(theta))) .* alt_alpha(theta, metric_handle);

end