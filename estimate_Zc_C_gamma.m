function params = estimate_Zc_C_gamma(SBP, DBP, MAP, HR, SV_ref, P1, ESP)
% Estimate Zc, C, gamma from SBP, DBP, MAP, HR and reference SV
% Now uses Ts derived from pressure morphology:
%   Ts = Ts_Weissler * (ESP - DBP) / (P1 - DBP)

% -----------------------------
% Prepare data
% -----------------------------
PP = SBP - DBP;

% Baseline Weissler systolic duration
Ts0 = 0.406 - 0.0017 .* HR;

% Morphology-based scaling factor
scale = (ESP - DBP) ./ (P1 - DBP);

% Individualised systolic duration
Ts = Ts0 .* scale;

% -----------------------------
% Vicorder-like model function
% -----------------------------
model = @(p,PP,Ts) 0.5 .* Ts .* ( PP./p(1) + p(2) .* (PP./(p(3).*Ts)) );
% p(1) = Zc
% p(2) = C
% p(3) = gamma

% -----------------------------
% Initial guesses
% -----------------------------
p0 = [0.12, 1.5, 0.35];

% -----------------------------
% Parameter bounds
% -----------------------------
lb = [0.01, 0.1, 0.05];
ub = [1.00, 5.0, 1.00];

% -----------------------------
% Fit parameters
% -----------------------------
opts = optimoptions('lsqcurvefit','Display','off');
p_est = lsqcurvefit(@(p,x) model(p,PP,Ts), p0, PP, SV_ref, lb, ub, opts);

% -----------------------------
% Compute fitted SV
% -----------------------------
SV_fit = model(p_est, PP, Ts);

% -----------------------------
% Regular R²
% -----------------------------
SSE = sum((SV_ref - SV_fit).^2);
SST = sum((SV_ref - mean(SV_ref)).^2);
R2 = 1 - SSE/SST;

% -----------------------------
% Mean difference & SD of difference
% -----------------------------
diffs = SV_fit - SV_ref;
mean_diff = mean(diffs);
sd_diff   = std(diffs);

% -----------------------------
% ICC(2,1)
% -----------------------------
X = [SV_ref(:), SV_fit(:)];
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

ICC_21 = (MS_subject - MS_error) / ...
         (MS_subject + (k - 1)*MS_error + k*(MS_rater - MS_error)/N);

% -----------------------------
% Observed vs Predicted Plot
% -----------------------------
figure; hold on; grid on;

scatter(SV_ref, SV_fit, 50, 'filled', 'MarkerFaceAlpha', 0.7);

p_reg = polyfit(SV_ref, SV_fit, 1);
xline_vals = linspace(min(SV_ref), max(SV_ref), 100);
plot(xline_vals, polyval(p_reg, xline_vals), 'r-', 'LineWidth', 2);

plot(xline_vals, xline_vals, 'k--', 'LineWidth', 1.5);

xlabel('Observed SV (reference)');
ylabel('Predicted SV (model)');
title('Observed vs Predicted Stroke Volume');

legend('Data', ...
       sprintf('Regression: y = %.2fx + %.2f', p_reg(1), p_reg(2)), ...
       'Identity line', ...
       'Location', 'best');

hold off;

% -----------------------------
% Output structure
% -----------------------------
params.Zc        = p_est(1);
params.C         = p_est(2);
params.gamma     = p_est(3);
params.R2        = R2;
params.mean_diff = mean_diff;
params.sd_diff   = sd_diff;
params.ICC_21    = ICC_21;
params.Ts        = Ts;

end