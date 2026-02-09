function stats = agree(x, y)
% compare_two_arrays
% Computes:
%   - mean difference (bias)
%   - standard deviation of the difference
%   - ICC(2,1) (two-way random, absolute agreement)
%
% INPUTS:
%   x, y : numeric vectors of equal length
%
% OUTPUT:
%   stats.mean_diff
%   stats.sd_diff
%   stats.icc

% Ensure column vectors
x = x(:);
y = y(:);

% -------------------------
% Mean difference & SD
% -------------------------
diffs = x - y;
stats.mean_diff = mean(diffs);
stats.sd_diff   = std(diffs);

% -------------------------
% ICC(2,1)
% -------------------------
X = [x y];   % N × 2 matrix
[N, k] = size(X);

mean_raters   = mean(X,1);
mean_subjects = mean(X,2);
grand_mean    = mean(X(:));

SS_total   = sum((X(:) - grand_mean).^2);
SS_subject = k * sum((mean_subjects - grand_mean).^2);
SS_rater   = N * sum((mean_raters - grand_mean).^2);
SS_error   = SS_total - SS_subject - SS_rater;

MS_subject = SS_subject / (N - 1);
MS_rater   = SS_rater   / (k - 1);
MS_error   = SS_error   / ((N - 1)*(k - 1));

stats.icc = (MS_subject - MS_error) / ...
            (MS_subject + (k - 1)*MS_error + ...
            k*(MS_rater - MS_error)/N);

end