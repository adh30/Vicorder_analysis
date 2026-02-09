function params = estimate_Zc_C_gamma_Ts(SBP, DBP, MAP, HR, SV_ref, P1, ESP)
% Fully integrated estimation of Zc, C, gamma, and Ts using:
% - 3-element Windkessel
% - Triangular inflow peaking at P1
% - ESP constraint to solve for Ts
% - Iterative refinement of Ts and (Zc, C, gamma)

N = length(SV_ref);

%% -----------------------------
% 1. Initial guesses
%% -----------------------------
Zc0    = 0.12;
C0     = 1.5;
gamma0 = 0.35;

% R = MAP / CO = MAP / (SV * HR/60)
CO = SV_ref .* HR/60;
R0 = median(MAP ./ CO);

alpha = 0.25;      % P1 occurs at 25% of Ts
P0    = DBP;       % pressure at foot

p = [Zc0, C0, gamma0, R0];

%% -----------------------------
% 2. Iterative loop
%% -----------------------------
for iter = 1:20

    Zc    = p(1);
    C     = p(2);
    gamma = p(3);
    R     = p(4);

    % --- Compute Ts for each beat ---
    Ts = zeros(N,1);
    for i = 1:N
        Ts(i) = solve_Ts_from_ESP(SV_ref(i), C, R, P0(i), alpha, ESP(i));
    end

    % Compute PP
    PP = SBP - DBP;

    % SV model
    SV_model = @(x) 0.5 .* Ts .* ( PP./x(1) + x(2) .* (PP./(x(3).*Ts)) );

    % Fit Zc, C, gamma
    fun = @(x) SV_model(x) - SV_ref;

    x0 = [Zc, C, gamma];
    lb = [0.01, 0.1, 0.05];
    ub = [1.00, 5.0, 1.00];

    opts = optimoptions('lsqnonlin','Display','off');
    x_est = lsqnonlin(fun, x0, lb, ub, opts);

    % Update parameters
    p(1:3) = x_est;

    % Update R
    CO = SV_ref .* HR/60;
    p(4) = median(MAP ./ CO);

end

%% -----------------------------
% 3. Final outputs
%% -----------------------------
Zc    = p(1);
C     = p(2);
gamma = p(3);
R     = p(4);

% Final Ts
Ts = zeros(N,1);
for i = 1:N
    Ts(i) = solve_Ts_from_ESP(SV_ref(i), C, R, P0(i), alpha, ESP(i));
end

PP = SBP - DBP;
SV_fit = 0.5 .* Ts .* ( PP./Zc + C .* (PP./(gamma.*Ts)) );

% Regular R2
SSE = sum((SV_ref - SV_fit).^2);
SST = sum((SV_ref - mean(SV_ref)).^2);
R2 = 1 - SSE/SST;

% Mean diff and SD
diffs = SV_fit - SV_ref;
mean_diff = mean(diffs);
sd_diff   = std(diffs);

% ICC(2,1)
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

% Output
params.Zc        = Zc;
params.C         = C;
params.gamma     = gamma;
params.R         = R;
params.Ts        = Ts;
params.SV_fit    = SV_fit;
params.R2        = R2;
params.mean_diff = mean_diff;
params.sd_diff   = sd_diff;
params.ICC_21    = ICC_21;

end