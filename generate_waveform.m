function [t, P, Qin] = generate_waveform(Zc, C, R, Ts, SV, P0, alpha)
% Generate a full pressure waveform using:
% - 3-element Windkessel (C, R)
% - Triangular inflow peaking at alpha*Ts
% - Stroke volume SV
% - Systolic duration Ts
%
% Inputs:
%   Zc    - characteristic impedance (not used in ODE, but included for completeness)
%   C     - compliance
%   R     - resistance
%   Ts    - systolic duration
%   SV    - stroke volume
%   P0    - initial pressure (DBP)
%   alpha - fractional timing of P1 (0 < alpha < 1)
%
% Outputs:
%   t     - time vector
%   P     - pressure waveform
%   Qin   - inflow waveform

%% Time discretisation
N  = 400;                 % number of steps
dt = Ts / N;
t  = linspace(0, Ts, N+1);

%% Triangular inflow peaking at t = alpha*Ts
Qpk = 2 * SV / Ts;        % from stroke volume constraint
Qin = zeros(size(t));

for k = 1:numel(t)
    tk = t(k);
    if tk <= alpha*Ts
        Qin(k) = Qpk * tk/(alpha*Ts);
    else
        Qin(k) = Qpk * (Ts - tk)/((1-alpha)*Ts);
    end
end

%% Forward Euler integration of Windkessel ODE
P = zeros(size(t));
P(1) = P0;

for k = 1:N
    dPdt = (Qin(k) - P(k)/R) / C;
    P(k+1) = P(k) + dt*dPdt;
end

end