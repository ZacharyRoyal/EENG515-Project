function samples = tri_samples(time_series, frequency)

    samples = sawtooth(frequency.*time_series, 0.5);
    
end