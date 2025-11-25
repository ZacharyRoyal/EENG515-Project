function samples = shark_samples(time_series, frequency)

    samples = tri_samples(time_series, frequency) + saw_samples(time_series, frequency);
    
end