function samples = saw_samples(time_series, frequency)

    samples = sawtooth(frequency.*time_series);
    
end