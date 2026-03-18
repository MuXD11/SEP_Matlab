%% BER vs Received Power — PIN photodiode, OOK modulation
clear; clc; close all;
 
% Physical constants
q   = 1.602e-19;        % [C]   electron charge
k_B = 1.381e-23;        % [J/K] Boltzmann constant
 
% Detector parameters
R   = 1;                % [A/W] responsivity
B   = 1e9;              % [Hz]  electrical bandwidth
T   = 290;              % [K]   noise temperature
F   = 3;                % [-]   noise figure (linear)
R_L = 50;               % [Ohm] load resistance
 
% Received power sweep
P_rx  = logspace(-12, -3, 500);        % [W]
P_dBm = 10 * log10(P_rx * 1e3);       % [dBm]
 
% Noise variances
var_shot = 2 * q * R .* P_rx * B;     % shot noise (signal-dependent)
var_th   = 4 * k_B * T * F * B / R_L; % thermal noise (constant)
 
% Q-factor and BER
Q   = R .* P_rx ./ (sqrt(var_shot + var_th) + sqrt(var_th));
BER = 0.5 * erfc(Q / sqrt(2));
 
% Plot
figure;
semilogy(P_dBm, BER, 'b-', 'LineWidth', 2);
yline(1e-9, 'r--', 'BER = 10^{-9}', 'LabelOrientation', 'horizontal');
xlabel('Received power (dBm)');
ylabel('BER');
title('BER vs Received Power');
grid on;
xlim([-100 -20]);
ylim([1e-12 1]);