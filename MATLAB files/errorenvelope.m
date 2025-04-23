% function to plot std error bars as envelopes
function errorenvelope(x, mean, stdev, color, FaceAlpha)
    % Ensure inputs are column vectors
    x = x(:);
    mean = mean(:);
    stdev = stdev(:);

    % Validate input dimensions
    assert(length(x) == length(mean) && length(mean) == length(stdev), ...
        'Inputs x, mean, and stdev must have the same length.');

    % Calculate upper and lower bounds
    y_low = mean - stdev;
    y_up = mean + stdev;

    % Prevent issues with log-log plots
    y_low = max(y_low, eps); % Avoid zero or negative values
    y_up = max(y_up, eps);

    % Create fill area for linear or log-log plots
    x_fill = [x; flipud(x)];
    y_fill = [y_up; flipud(y_low)];

    % Plot the envelope
    hFill = fill(x_fill, y_fill, color, 'FaceAlpha', FaceAlpha, 'LineStyle', 'none', 'EdgeColor', 'none');

    % Optionally remove from legend
    set(get(get(hFill, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end