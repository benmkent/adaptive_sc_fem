function [x_mean, x_std, x_min, x_max] = extract_summary_stats(x)
    x_mean = mean(x);
    x_std = std(x);
    x_min = min(x);
    x_max = max(x);
end