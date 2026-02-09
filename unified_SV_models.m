function SV = unified_SV_models(SBP, DBP, MAP, HR)
% unified_SV_models_pressure_only
% Computes stroke volume using four toy model families:
%   1. Modelflow-like (pressure-only)
%   2. Vicorder-like (pressure-only)
%   3. ARCSolver-like (pressure-only)
%   4. Liljestran-Zander
%   5. PiCCO
%   6. PP

% -----------------------------
% Shared quantities
% -----------------------------
PP = SBP - DBP;          % pulse pressure
T  = 60 ./ HR;           % cardiac cycle length (s)
% Ts = 0.3 .* T;         % systolic duration (approx)
Ts = 0.33 + 0.0012.*HR;  % Weissler's empirical formula
Td = T-Ts;               % diastolic duration 

% Systolic dP/dt estimate (scalar)
dPdt = 1.5*(SBP - DBP) ./ Ts;

% ============================================================
% 1. MODELFLOW-LIKE (BEST ESTIMATE)
% ============================================================
% Zc_mf = 0.12;                           % mmHg·s/mL
Zc_mf = 0.12;                           % aortic values 0.15?
C_mf  = .6465 .* exp(-0.02 .* (MAP - 90)); % 1, 0.02 and 90 are approximate constants

Qpk_mf = PP ./ Zc_mf + C_mf .* dPdt;
SV.modelflow = 0.5 .* Qpk_mf .* Ts;

% ============================================================
% 2. VICORDER-LIKE (ORIGINAL)
% ============================================================
Zc_v = 0.14;
C_v  = 1.5 .* exp(-0.025 .* (MAP - 80));

Qpk_v = PP ./ Zc_v + C_v .* dPdt;
SV.vicorder = 0.5 .* Qpk_v .* Ts;

% ============================================================
% 3. ARCSOLVER-LIKE 
% ============================================================
C_a = 1 .* exp(-0.025 .* (MAP - 80));
R_a = MAP ./ (5.0 .* 1000 ./ 60);       % mmHg·s/mL

Qpk_a = C_a .* dPdt + PP ./ R_a;
SV.arcsolver = .5 .* Qpk_a .* Ts;

% ============================================================
% 4. Liljestran-Zander
% ============================================================
k_lz = 1000/1.331;
SV.lz = k_lz*PP./(SBP+DBP);

% ============================================================
% 5. PiCCO
% ============================================================
k_picco = 6.4103;    % empirical
sum_sbp=Ts.*SBP*.7; % approximation
SV.picco = k_picco*(1+(Ts./Td).*sum_sbp);

% ============================================================
% 6. PP Improved Pulse Pressure Model
% ============================================================
Zc_pp = 0.12;
C_pp = 0.6465;
gamma = 0.7476;
triang = 0.5;  
SV.pp = triang*Ts.*(PP./Zc_pp+C_pp.*(PP./(gamma.*Ts)));

end