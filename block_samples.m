function samples = block_samples(time_series, frequency)

    samples = floor(sin(frequency.*time_series))+0.5;
    
end