function [D, dict_info] = build_alt_dictionary(N, freqs, metric_def_list)
    % Concatenates columns from all metrics/frequencies
    cols = {};
    dict_info = {};
    for m = 1:numel(metric_def_list)
        Phi_m = build_alt_basis(N, freqs, metric_def_list{m});
        for j = 1:numel(freqs)
            cols{end+1} = Phi_m(:, j);
            dict_info{end+1} = struct('freq', freqs(j), 'metric', metric_def_list{m});
        end
    end
    D = cell2mat(cols); % N x K
end
