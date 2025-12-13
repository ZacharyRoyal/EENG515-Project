function [D, dict_info] = build_alt_dictionary(N, metric_def_list)
    % Concatenates columns from all metrics/frequencies
    freqs = 0:N-1;
    cols = {};
    dict_info = {};
    for m = 1:numel(metric_def_list)
        Phi_m = build_alt_basis(N, metric_def_list{m});
        for j = 1:numel(freqs)
            cols{end+1} = Phi_m(:, j);
            dict_info{end+1} = struct('freq', freqs(j), 'metric', metric_def_list{m});
        end
    end
    D = cell2mat(cols); % N x K
    fprintf("Size of Dictionary, K: %d\n", size(D, 2));
end
