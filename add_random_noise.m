function samples = add_random_noise(samples, noise_power_scalar, smoothing_level)

    noise = smoothed_random_samples(length(samples),smoothing_level).*noise_power_scalar;

    samples = samples + noise;
    
end