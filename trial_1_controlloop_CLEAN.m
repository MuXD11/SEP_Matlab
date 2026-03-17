%% FULL POINTING-JITTER / NORMALIZED-POWER SIMULATION
% This script implements the complete chain:
% 1) FSM plant P(s) = 1 / (m s^2 + c s + k)
% 2) PID-like controller C(s)
% 3) Random disturbance generation in x and y
% 4) Residual disturbance through sensitivity S = 1/(1+PC)
% 5) Gaussian beam mapping to normalized power I(t)
%
% It produces:
% - open-loop disturbance signals x_d(t), y_d(t)
% - residual closed-loop disturbance x_r(t), y_r(t)
% - normalized received power traces I_open(t), I_closed(t)

clear; clc; close all;

% ------------------------------------------------------------------------
% 1) PARAMETERS
% -------------------------------------------------------------------------

% Mechanical FSM / mass-spring-damper parameters
m = 0.1;        % kg
c = 1.0;        % Ns/m
k = 400.0;      % N/m

% Controller design target
fc = 500;               % target crossover frequency [Hz]
wc = 2*pi*fc;           % rad/s

% Controller shaping parameters
wd  = 9;                % derivative-shaping factor
wi  = 10;               % integral-shaping factor
wis = 1;                % leaky integrator pole factor

% Disturbance parameters
f_bw = 100;             % disturbance bandwidth [Hz]
sigma_dist = 10e-6;     % target std of disturbance [rad] = 10 urad

% Optical beam parameters
w_beam = 1e-4;          % beam width parameter [rad]
x_se = 0;               % static x bias / boresight error [rad]
y_se = 0;               % static y bias / boresight error [rad]

% TODO: Play with systematic error
% TODO: Smaller angular divergence

% Simulation settings
fs = 10*fc;             % simulation sample frequency [Hz]
dt = 1/fs;
Tsim = 5.0;             % seconds
t = (0:dt:Tsim).';      % column vector time axis
N = numel(t);

rng(1);                 % reproducible random signals

% ------------------------------------------------------------------------
% 2) BUILD THE FSM PLANT P(s)
% -------------------------------------------------------------------------
% Plant: P(s) = x(s)/F(s) = 1 / (m s^2 + c s + k)

s = tf('s');
P_fsm = 1 / (m*s^2 + c*s + k);

% ------------------------------------------------------------------------
% 3) BUILD THE CONTROLLER C(s)
% -------------------------------------------------------------------------
% C(s) = Kp + Ki/(s + 1/taui) + Kd*s/(Kt*s + 1)

Kp   = m*wc^2 / sqrt(wd);
Kd   = Kp*sqrt(wd) / wc;
Kt   = 1 / (sqrt(wd)*wc);
Ki   = Kp*(wc/wi);
taui = 1 / (wis*2*pi); % Integral filterting

C = (Kp + Ki*tf(1,[1 1/taui])) + Kd*tf([1 0],[Kt 1]);

% ------------------------------------------------------------------------
% 4) OPEN-LOOP AND SENSITIVITY. minimal realization finds & removes cancelling pole-zero pairs
% -------------------------------------------------------------------------
L = minreal(C * P_fsm);       % open-loop transfer function L = P*C
S = minreal(1 / (1 + L));     % sensitivity:               S = 1/(1+L)
T = minreal(L / (1 + L));     % complementary sensitivity: T = L/(1+L)

% ------------------------------------------------------------------------
% 5) GENERATE RANDOM DISTURBANCE SOURCES
% -------------------------------------------------------------------------
nx = randn(N,1);
ny = randn(N,1);

% ------------------------------------------------------------------------
% 6) SHAPE THE DISTURBANCES WITH A 2nd-ORDER FILTER
% -------------------------------------------------------------------------
[numF, denF] = butter(2, 2*pi*f_bw, 'low', 's');
P_noise = tf(numF, denF);

xd_raw = lsim(P_noise, nx, t);% Convolution
yd_raw = lsim(P_noise, ny, t);

% Normalize to the desired standard deviation = 10 urad (sigma_dist). Standar dev scales
% linearly 
xd = xd_raw * sigma_dist / std(xd_raw);
yd = yd_raw * sigma_dist / std(yd_raw);

% ------------------------------------------------------------------------
% 7) APPLY CLOSED-LOOP DISTURBANCE REJECTION
% -------------------------------------------------------------------------
xr = lsim(S, xd, t);
yr = lsim(S, yd, t);

% ------------------------------------------------------------------------
% 8) COMPUTE NORMALIZED RECEIVED POWER
% -------------------------------------------------------------------------
I_open   = exp(-2 * ((xd + x_se).^2 + (yd + y_se).^2) / w_beam^2);
I_closed = exp(-2 * ((xr + x_se).^2 + (yr + y_se).^2) / w_beam^2);

% ------------------------------------------------------------------------
% 9) PRINT USEFUL STATISTICS
% -------------------------------------------------------------------------
fprintf('--- Disturbance statistics ---\n');
fprintf('std(xd)  = %.3f urad\n', std(xd)*1e6);
fprintf('std(yd)  = %.3f urad\n', std(yd)*1e6);
fprintf('std(xr)  = %.3f urad\n', std(xr)*1e6);
fprintf('std(yr)  = %.3f urad\n', std(yr)*1e6);

fprintf('\n--- Power statistics ---\n');
fprintf('mean(I_open)   = %.4f\n', mean(I_open));
fprintf('mean(I_closed) = %.4f\n', mean(I_closed));
fprintf('min(I_open)    = %.4f\n', min(I_open));
fprintf('min(I_closed)  = %.4f\n', min(I_closed));

% ========================================================================
% 10) PLOTS
% ========================================================================

% Shared frequency vector for all Bode-style plots
% Range: 0.1 Hz to 10x crossover frequency
f_vec = logspace(-1, log10(10*fc), 1000);   % Hz
w_vec = 2*pi*f_vec;                          % rad/s

% Helper: evaluate a tf and return magnitude in dB and phase in degrees
%   mag_db, phase_deg are column vectors matching f_vec
eval_tf = @(G) deal( ...
    20*log10(squeeze(abs(freqresp(G, w_vec)))), ...
    rad2deg(squeeze(angle(freqresp(G, w_vec)))) );

% -----------------------------------------------------------------------
% PLOT A: FSM Plant P(s)
% -----------------------------------------------------------------------
[P_mag, P_pha] = eval_tf(P_fsm);

figure('Name','Plant P(s)');
subplot(2,1,1);
semilogx(f_vec, P_mag, 'b', 'LineWidth', 1.5);
grid on;
ylabel('Magnitude (dB)');
title('FSM Plant  P(s) = 1/(ms^2 + cs + k)');
yline(0, 'k--', '0 dB', 'LabelHorizontalAlignment','left');

subplot(2,1,2);
semilogx(f_vec, P_pha, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
ylim([-200 10]);
yticks([-180 -135 -90 -45 0]);

% -----------------------------------------------------------------------
% PLOT B: Open-Loop L(s) = P(s)*C(s)
%   Key property: |L(fc)| = 1  <=>  0 dB at crossover frequency fc
% -----------------------------------------------------------------------
[L_mag, L_pha] = eval_tf(L);

% Verify crossover numerically
[~, PM, ~, Wcp] = margin(L);
fc_actual = Wcp / (2*pi);
fprintf('\n--- Loop verification ---\n');
fprintf('Designed crossover: fc = %d Hz\n', fc);
fprintf('Actual  crossover: fc = %.1f Hz  (where |L|=1, i.e. 0 dB)\n', fc_actual);
fprintf('Phase margin: %.1f deg\n', PM);

figure('Name','Open-Loop L(s) = P(s)*C(s)');
subplot(2,1,1);
semilogx(f_vec, L_mag, 'b', 'LineWidth', 1.5);
hold on;
% Mark the 0 dB line — this is where the crossover happens
yline(0, 'k--', 'LineWidth', 1.2);
% Mark the actual crossover frequency with a vertical line
xline(fc_actual, 'r--', sprintf('f_c = %.0f Hz', fc_actual), ...
      'LineWidth', 1.2, 'LabelVerticalAlignment','bottom');
grid on;
ylabel('Magnitude (dB)');
title('Open-Loop  L(s) = P(s) \cdot C(s)');
legend('|L(j\omega)|', '0 dB (crossover line)', ...
       sprintf('Crossover f_c = %.0f Hz', fc_actual), ...
       'Location','southwest');

subplot(2,1,2);
semilogx(f_vec, L_pha, 'b', 'LineWidth', 1.5);
hold on;
yline(-180, 'k--', '-180 deg', 'LabelHorizontalAlignment','left');
xline(fc_actual, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
ylim([-270 10]);
yticks([-270 -180 -135 -90 -45 0]);

% -----------------------------------------------------------------------
% PLOT C: Sensitivity S(s) and Complementary Sensitivity T(s)
%   S + T = 1  (always true — they are complementary)
%   S small → disturbances rejected  (inside bandwidth)
%   T large → reference tracked well (inside bandwidth)
% -----------------------------------------------------------------------
[S_mag, S_pha] = eval_tf(S);
[T_mag, T_pha] = eval_tf(T);

figure('Name','S(s) and T(s)');
subplot(2,1,1);
semilogx(f_vec, S_mag, 'b',  'LineWidth', 1.5); hold on;
semilogx(f_vec, T_mag, 'r',  'LineWidth', 1.5);
yline(0,  'k--', '0 dB',  'LabelHorizontalAlignment','left');
xline(fc_actual, 'k:', sprintf('f_c = %.0f Hz', fc_actual), ...
      'LineWidth', 1.0, 'LabelVerticalAlignment','bottom');
grid on;
ylabel('Magnitude (dB)');
title('Sensitivity S(s) and Complementary Sensitivity T(s)');
legend('S = 1/(1+L)  [disturbance rejection]', ...
       'T = L/(1+L)  [reference tracking]', ...
       '0 dB', 'Crossover f_c', ...
       'Location','east');
% Annotation: S small = good rejection, T~1 = good tracking
text(0.15, 0.15, ...
    {'S \approx 0: disturbance rejected'; 'T \approx 1: reference tracked'}, ...
    'Units','normalized','FontSize',8,'Color',[0 0.5 0]);
text(0.65, 0.45, ...
    {'S \approx 1: disturbance'; 'passes through'}, ...
    'Units','normalized','FontSize',8,'Color',[0.7 0 0]);

subplot(2,1,2);
semilogx(f_vec, S_pha, 'b', 'LineWidth', 1.5); hold on;
semilogx(f_vec, T_pha, 'r', 'LineWidth', 1.5);
xline(fc_actual, 'k:', 'LineWidth', 1.0);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
legend('Phase of S', 'Phase of T', 'Location','best');

% -----------------------------------------------------------------------
% PLOT D: Disturbances before and after control
% -----------------------------------------------------------------------
figure('Name','Disturbances: open vs closed loop');
subplot(2,1,1);
plot(t, xd*1e6, 'r', 'LineWidth', 0.8); hold on;
plot(t, xr*1e6, 'b', 'LineWidth', 0.8);
grid on;
ylabel('x (\murad)');
title('X disturbance: open-loop vs closed-loop');
legend(sprintf('Open-loop  std = %.2f \\murad', std(xd)*1e6), ...
       sprintf('Closed-loop std = %.3f \\murad', std(xr)*1e6));

subplot(2,1,2);
plot(t, yd*1e6, 'r', 'LineWidth', 0.8); hold on;
plot(t, yr*1e6, 'b', 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('y (\murad)');
title('Y disturbance: open-loop vs closed-loop');
legend('Open-loop', 'Closed-loop');

% -----------------------------------------------------------------------
% PLOT E: Normalized received power
% -----------------------------------------------------------------------
figure('Name','Normalized received power');
plot(t, I_open,   'r', 'LineWidth', 0.8); hold on;
plot(t, I_closed, 'b', 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Normalized power  I(t)');
title('Normalized received power: open-loop vs closed-loop');
legend(sprintf('Open-loop   mean = %.3f', mean(I_open)), ...
       sprintf('Closed-loop mean = %.4f', mean(I_closed)), ...
       'Location','best');
ylim([0 1.05]);

% -----------------------------------------------------------------------
% PLOT F: Power distributions
% -----------------------------------------------------------------------
figure('Name','Distribution of normalized power');
histogram(I_open,   80, 'Normalization','pdf', ...
          'FaceColor','r', 'FaceAlpha',0.5); hold on;
histogram(I_closed, 80, 'Normalization','pdf', ...
          'FaceColor','b', 'FaceAlpha',0.5);
grid on;
xlabel('Normalized power  I');
ylabel('PDF estimate');
title('Distribution of normalized power');
legend('Open-loop', 'Closed-loop');