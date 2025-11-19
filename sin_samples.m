function samples = sin_samples(time_series, frequency, metric_handle)

    samples = alt_sin(frequency.*time_series, metric_handle);
    
end