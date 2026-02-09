function [Ts, diagTs] = solve_Ts_from_pressures(SV, C, R, P0, alpha, ...
                                                P1_obs, SBP_obs, ESP_obs, ...
                                                doPlot)
% Estimate Ts by minimising error between modelled and observed:
%   - P1
%   - SBP
%   - ESP
% using a 3-element Windkessel + triangular inflow peaking at P1.

if nargin < 10
    doPlot = false;
end

tau = R * C;

% Physiological bounds
Ts_lower = 0.20;
Ts_upper = 0.45;

% Weighting of pressure constraints
wP1  = 1.0;
wSBP = 1.0;
wESP = 1.0;

% Objective function
obj = @(Ts) objective_Ts(Ts, SV, C, R, P0, alpha, ...
                         P1_obs, SBP_obs, ESP_obs, ...
                         wP1, wSBP, wESP, tau);

% Minimise over Ts
opts = optimset('Display','off');
[Ts_opt, Jmin] = fminbnd(obj, Ts_lower, Ts_upper, opts);

Ts = Ts_opt;

% Diagnostics: boundary-limited?
epsB = 0.005;
boundary_flag = (abs(Ts - Ts_lower) < epsB) | (abs(Ts - Ts_upper) < epsB);

diagTs.boundary_flag = boundary_flag;
diagTs.Jmin          = Jmin;

% Optional diagnostic plot
if doPlot
    Ts_grid = linspace(Ts_lower, Ts_upper, 200);
    ESP_grid = zeros(size(Ts_grid));
    for k = 1:numel(Ts_grid)
        [~, ~, ESP_grid(k)] = forward_pressures(Ts_grid(k), SV, C, R, P0, alpha, tau);
    end

    figure; hold on;
    plot(Ts_grid, ESP_grid, 'b-', 'LineWidth', 1.5);
    yline(ESP_obs, 'r--', 'LineWidth', 1.5);
    xline(Ts, 'k:', 'LineWidth', 1.5);
    xlabel('T_s (s)');
    ylabel('ESP_{model} (mmHg)');
    legend('ESP_{model}(T_s)', 'ESP_{obs}', 'T_s^*');
    title('ESP constraint vs T_s');
    grid on;
end

end


%% ============================================================
%  Objective function for Ts
%% ============================================================
function J = objective_Ts(Ts, SV, C, R, P0, alpha, ...
                          P1_obs, SBP_obs, ESP_obs, ...
                          wP1, wSBP, wESP, tau)

[P1_mod, SBP_mod, ESP_mod] = forward_pressures(Ts, SV, C, R, P0, alpha, tau);

eP1  = P1_mod  - P1_obs;
eSBP = SBP_mod - SBP_obs;
eESP = ESP_mod - ESP_obs;

J = wP1*eP1.^2 + wSBP*eSBP.^2 + wESP*eESP.^2;

end


%% ============================================================
%  Forward Windkessel model (numerical)
%% ============================================================
function [P1_mod, SBP_mod, ESP_mod] = forward_pressures(Ts, SV, C, R, P0, alpha, tau)

% Time discretisation
N  = 400;
dt = Ts / N;
t  = linspace(0, Ts, N+1);

% Triangular inflow peaking at t = alpha*Ts
Qpk = 2 * SV / Ts;
Qin = zeros(size(t));

for k = 1:numel(t)
    tk = t(k);
    if tk <= alpha*Ts
        Qin(k) = Qpk * tk/(alpha*Ts);
    else
        Qin(k) = Qpk * (Ts - tk)/((1-alpha)*Ts);
    end
end

% Forward Euler integration of Windkessel ODE
P = zeros(size(t));
P(1) = P0;

for k = 1:N
    dPdt = (Qin(k) - P(k)/R) / C;
    P(k+1) = P(k) + dt*dPdt;
end

% Extract model pressures
[~, idxP1] = min(abs(t - alpha*Ts));
P1_mod  = P(idxP1);
SBP_mod = max(P);
ESP_mod = P(end);

end