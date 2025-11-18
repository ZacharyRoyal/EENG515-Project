% takes in raw dft coefficients and plots them nicely
function plot_disc_fourier(coeffs)

    % extract the magnitudes
    mags = abs(coeffs);

    % surpress teeny tiny coefficients
    surpressed_coeffs = coeffs;
    surpressed_coeffs(mags<1e-6) = 0;

    % and then extract the phases 
    phases = unwrap(angle(surpressed_coeffs));

    % create vector of frequencies to plot over
    freqs = ((0:length(coeffs)-1)*100)/length(coeffs);

    subplot(2,1,1)
    plot(freqs, mags)
    title('Magnitude')

    subplot(2,1,2)
    plot(freqs, phases*180/pi)
    title('Phase')

end