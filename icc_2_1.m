function ICC = icc_2_1(X)
% ICC(2,1): Two-way random effects, absolute agreement
% X is an N×k matrix (N subjects, k raters/methods)

[N, k] = size(X);

% Means
mean_raters = mean(X,1);
mean_subjects = mean(X,2);
grand_mean = mean(X(:));

% Sum of squares
SS_total = sum((X(:) - grand_mean).^2);
SS_between_subjects = k * sum((mean_subjects - grand_mean).^2);
SS_between_raters = N * sum((mean_raters - grand_mean).^2);
SS_residual = SS_total - SS_between_subjects - SS_between_raters;

% Mean squares
MS_subjects = SS_between_subjects / (N - 1);
MS_raters   = SS_between_raters / (k - 1);
MS_error    = SS_residual / ((N - 1)*(k - 1));

% ICC(2,1)
ICC = (MS_subjects - MS_error) / ...
      (MS_subjects + (k - 1)*MS_error + k*(MS_raters - MS_error)/N);
end