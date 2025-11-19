function samples = smoothed_random_samples(n_samples, smoothing_level)

    samples = randn(1,n_samples);
    windowSize = max([round(n_samples/20),3]); 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;

    for i = 1:1:smoothing_level
        samples = filter(b,a,samples);
    end
    
end