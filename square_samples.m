function samples = square_samples(time_series, frequency)

    samples = round(sin(frequency.*time_series));
    
end