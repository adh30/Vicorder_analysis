# Vicorder_analysis
# Windkessel-Based Stroke Volume Estimation

This repository contains MATLAB code for beat-by-beat estimation of
arterial parameters using a three-element Windkessel model with a
triangular inflow waveform peaking at the first systolic shoulder (P1).

## Features

- Beat-by-beat estimation of:
  - Systolic duration (Ts)
  - Characteristic impedance (Zc)
  - Arterial compliance (C)
  - Shape parameter (gamma)
  - Total peripheral resistance (R)
  - Model-predicted stroke volume (SV_fit)
- Numerical Windkessel forward model
- Triangular inflow waveform constrained by stroke volume
- Regularised nonlinear fitting of Zc, C, gamma
- Export of all parameters to Excel
- Optional waveform generation

## File Overview

### `estimate_Zc_C_gamma_Ts_perBeat.m`
Main driver function. Computes Ts, Zc, C, gamma, R, and SV_fit for each beat.

### `solve_Ts_from_pressures.m`
Estimates systolic duration Ts by minimising error between modelled and
observed P1, SBP, and ESP.

### `forward_pressures.m`
Numerical Windkessel forward model used inside the Ts solver.

### `generate_waveform.m`
Generates a full synthetic pressure waveform for a given beat.

### `export_params_to_excel.m`
Exports the parameter structure to an Excel file.

## Usage Example

```matlab
params = estimate_Zc_C_gamma_Ts_perBeat( ...
    SBP, DBP, MAP, HR, SV_ref, P1, ESP);

export_params_to_excel(params, 'Windkessel_Params.xlsx');
