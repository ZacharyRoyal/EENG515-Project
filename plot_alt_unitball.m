function plot_alt_unitball(metric_handle)
    % Number of points
    N = 10000;
    theta = linspace(0, 2*pi, N);

    % Parametrize the boundary of the p-norm ball
    x = cos(theta);
    y = sin(theta);

    % Normalize each point to lie on the p-norm unit ball
    alpha = 1 ./ metric_handle(x, y);
    x = x .* alpha;
    y = y .* alpha;

    % Plot the unit ball boundary
    figure; 
    subplot(1,3,1);
    plot(x, y, 'b-', 'LineWidth', 2);
    axis equal;
    title(['Unit ball under metric:', func2str(metric_handle)]);
    xlabel('x'); ylabel('y');

    % Plot "generalized sine" and "generalized cosine"
    subplot(1,3,2);
    plot(theta, x, 'r-', 'LineWidth', 1.5); hold on;
    plot(theta, y, 'g-', 'LineWidth', 1.5);
    legend('cos_{alt}(\theta)', 'sin_{alt}(\theta)');
    title('Generalized sine and cosine vs. angle parameter');
    xlabel('\theta'); ylabel('value');

    %Plot alpha
    subplot(1,3,3);
    plot(theta, alpha, 'b-', 'LineWidth', 1.5); hold on;
    legend('\alpha(\theta)');
    title('\alpha(\theta), scalar applied to base sin and cos');
    xlabel('\theta'); ylabel('value');

end

