function Ts = solve_Ts_from_ESP(SV, C, R, P0, alpha, ESP_obs)

tau = R * C;

ESP_model = @(Ts) compute_ESP(Ts, SV, C, R, P0, alpha, tau);
f = @(Ts) ESP_model(Ts) - ESP_obs;

Ts_lower = 0.20;
Ts_upper = 0.45;

fL = f(Ts_lower);
fU = f(Ts_upper);

% --- Case 1: Root exists in interval ---
if fL * fU < 0
    Ts = fzero(f, [Ts_lower, Ts_upper]);
    return
end

% --- Case 2: No root in interval ---
% Choose Ts that minimises |f(Ts)|

Ts_candidates = linspace(Ts_lower, Ts_upper, 200);
[~, idx] = min(abs(f(Ts_candidates)));
Ts = Ts_candidates(idx);

end


function ESP = compute_ESP(Ts, SV, C, R, P0, alpha, tau)

% I2
a = (1 - alpha) * Ts;
I2 = (tau^2 - (tau*a + tau^2).*exp(-a./tau)) ./ ((1 - alpha)*Ts);

% I1
term1 = (tau*alpha*Ts - tau^2).*exp(-(1 - alpha)*Ts./tau);
term2 = tau^2.*exp(-Ts./tau);
I1 = (term1 + term2) ./ (alpha*Ts);

ESP = P0.*exp(-Ts./tau) + (2*SV./(C*Ts)).*(I1 + I2);

end