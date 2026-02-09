function params = estimate_Zc_C_gamma_Ts_perBeat(SBP, DBP, MAP, HR, SV_ref, P1, ESP)

N = length(SV_ref);

% Outputs
Zc      = zeros(N,1);
C       = zeros(N,1);
gamma   = zeros(N,1);
R       = zeros(N,1);
Ts      = zeros(N,1);
SV_fit  = zeros(N,1);
Ts_flag = false(N,1);

% Regularisation priors
Zc0    = 0.12;
C0     = 1.5;
gamma0 = 0.35;
lambda = 0.1;

alpha = 0.25;   % P1 occurs at 25% of Ts

for i = 1:N

    %% -----------------------------
    % 1. Estimate R for this beat
    %% -----------------------------
    CO_i = SV_ref(i) * HR(i) / 60;   % cardiac output for beat i
    R_i  = MAP(i) / CO_i;            % resistance for beat i
    R(i) = R_i;

    %% -----------------------------
    % 2. Estimate Ts for this beat
    %% -----------------------------
    P0_i = DBP(i);

    [Ts_i, diagTs] = solve_Ts_from_pressures( ...
        SV_ref(i), C0, R_i, P0_i, alpha, ...
        P1(i), SBP(i), ESP(i), false);

    Ts(i)      = Ts_i;
    Ts_flag(i) = diagTs.boundary_flag;

    %% -----------------------------
    % 3. Fit Zc, C, gamma for this beat
    %% -----------------------------
    PP_i = SBP(i) - DBP(i);

    SV_model = @(x) 0.5 .* Ts_i .* ( PP_i./x(1) + x(2) .* (PP_i./(x(3).*Ts_i)) );

    fun = @(x) [ SV_model(x) - SV_ref(i); ...
                 sqrt(lambda)*(x(1) - Zc0); ...
                 sqrt(lambda)*(x(2) - C0);  ...
                 sqrt(lambda)*(x(3) - gamma0) ];

    x0 = [Zc0, C0, gamma0];
    lb = [0.01, 0.1, 0.05];
    ub = [1.00, 5.0, 1.00];

    opts = optimoptions('lsqnonlin','Display','off');
    x_est = lsqnonlin(fun, x0, lb, ub, opts);

    Zc(i)    = x_est(1);
    C(i)     = x_est(2);
    gamma(i) = x_est(3);

    %% -----------------------------
    % 4. Compute fitted SV
    %% -----------------------------
    SV_fit(i) = SV_model(x_est);

end

%% -----------------------------
% 5. Compute global diagnostics
%% -----------------------------
SSE = sum((SV_ref - SV_fit).^2);
SST = sum((SV_ref - mean(SV_ref)).^2);
R2  = 1 - SSE/SST;

diffs     = SV_fit - SV_ref;
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

%% -----------------------------
% 6. Output
%% -----------------------------
params.Zc        = Zc;
params.C         = C;
params.gamma     = gamma;
params.R         = R;
params.Ts        = Ts;
params.SV_fit    = SV_fit;
params.Ts_flag   = Ts_flag;

params.R2        = R2;
params.mean_diff = mean_diff;
params.sd_diff   = sd_diff;
params.ICC_21    = ICC_21;

end