%% BER vs SNR — Physical detector noise + pointing jitter fading
% Merges Alice's PIN photodiode noise model with Carlos's fading simulation.
%
% Depends on: trial_1_controlloop_CLEAN.m  (run first to get I_closed)
%
% Outer loop sweeps average received optical power P0.
% For each P0, shot noise varies per fading block, thermal noise is fixed.
% x-axis is the SNR derived from the physics, not a free parameter.

% -------------------------------------------------------------------------
% 0) SIMULATION SETTINGS  (must match control-loop script)
% -------------------------------------------------------------------------
fs   = 5000;                % [Hz]  sample frequency
Tsim = 10;                  % [s]   BER simulation duration
N    = Tsim * fs;           % number of fading blocks

bit_rate        = 1e8;              % [bps] 100 Mbps
bits_per_sample = round(bit_rate/fs); % bits per fading block (= 20000)

% Trim I_closed to BER simulation length
I_sim = I_closed(1:N);

% -------------------------------------------------------------------------
% 1) DETECTOR PARAMETERS  (from Alice's code)
% -------------------------------------------------------------------------
q   = 1.602e-19;        % [C]    electron charge
k_B = 1.381e-23;        % [J/K]  Boltzmann constant
R   = 1;                % [A/W]  responsivity
B   = 1e8;              % [Hz]   electrical bandwidth = bit rate
T   = 290;              % [K]    noise temperature
F   = 3;                % [-]    noise figure (linear)
R_L = 100;              % [Ohm]  load resistance

% Thermal noise variance — fixed, independent of signal
var_th = 4 * k_B * T * F * B / R_L;   % [A^2]
sig_th = sqrt(var_th);                 % [A]

% -------------------------------------------------------------------------
% 2) RECEIVED POWER SWEEP
% -------------------------------------------------------------------------
% Sweep P0 over a physically meaningful range.
% At R=1 A/W and these noise levels, sensitivity is typically around -40 dBm.
P0_array  = logspace(-9, -3, 40);          % [W]  1 nW to 1 mW
P0_dBm    = 10*log10(P0_array * 1e3);      % [dBm]
num_P0    = numel(P0_array);

% -------------------------------------------------------------------------
% 3) PREALLOCATE RESULTS
% -------------------------------------------------------------------------
BER_MC      = zeros(1, num_P0);
BER_new     = zeros(1, num_P0);   % analytical: averaged over fading
BER_old     = zeros(1, num_P0);   % analytical: mean-I only (Alice's Q at mean)
SNR_dB_axis = zeros(1, num_P0);   % SNR derived from physics, for x-axis

% -------------------------------------------------------------------------
% 4) MAIN LOOP OVER P0 LEVELS
% -------------------------------------------------------------------------
for p_idx = 1:num_P0

    P0 = P0_array(p_idx);

    % ------------------------------------------------------------------
    % 4a) Per-block physical quantities
    % ------------------------------------------------------------------
    % Instantaneous received power for each fading block
    P_k = P0 .* I_sim;                        % [W]  Nx1 vector

    % Photocurrent for bit=1 at each block
    i1_k = R .* P_k;                          % [A]

    % Shot noise variance for bit=1 (signal-dependent)
    var_shot_k = 2 * q * R .* P_k * B;        % [A^2]
    sig1_k     = sqrt(var_shot_k + var_th);   % [A]  total noise std, bit=1

    % Bit=0: no signal, only thermal noise
    sig0 = sig_th;                             % [A]  scalar

    % ------------------------------------------------------------------
    % 4b) Optimal threshold per block  [A]
    %   tau(k) = sigma0 * i1(k) / (sigma1(k) + sigma0)
    %   Derived from crossing point of the two conditional Gaussians.
    % ------------------------------------------------------------------
    tau_k = (sig0 .* i1_k) ./ (sig1_k + sig0);   % [A]  Nx1 vector

    % ------------------------------------------------------------------
    % 4c) Analytical BER — averaged over fading (BER_new)
    %   Q(k) = i1(k) / (sigma1(k) + sigma0)
    %   Pb(k) = 0.5 * erfc(Q(k) / sqrt(2))
    % ------------------------------------------------------------------
    Q_k        = i1_k ./ (sig1_k + sig0);
    Pb_k       = 0.5 * erfc(Q_k / sqrt(2));
    BER_new(p_idx) = mean(Pb_k);

    % ------------------------------------------------------------------
    % 4d) Analytical BER — mean-I only (BER_old, Alice's formula at mean P)
    % ------------------------------------------------------------------
    P_mean     = P0 * mean(I_sim);
    i1_mean    = R * P_mean;
    sig1_mean  = sqrt(2*q*R*P_mean*B + var_th);
    Q_mean     = i1_mean / (sig1_mean + sig0);
    BER_old(p_idx) = 0.5 * erfc(Q_mean / sqrt(2));

    % ------------------------------------------------------------------
    % 4e) SNR at mean received power — used as x-axis
    %   SNR = i1_mean^2 / (var_shot_mean + var_th)
    % ------------------------------------------------------------------
    var_shot_mean     = 2 * q * R * P_mean * B;
    SNR_lin           = i1_mean^2 / (var_shot_mean + var_th);
    SNR_dB_axis(p_idx) = 10*log10(SNR_lin);

    % ------------------------------------------------------------------
    % 4f) Monte Carlo BER
    % ------------------------------------------------------------------
    total_errors = 0;
    total_bits   = 0;

    fprintf('[P0 = %.1f dBm | SNR = %.1f dB]\n', P0_dBm(p_idx), SNR_dB_axis(p_idx));

    for k = 1:N

        i1   = i1_k(k);      % signal current this block  [A]
        s1   = sig1_k(k);    % noise std for bit=1        [A]
        tau  = tau_k(k);     % optimal threshold           [A]

        for bit_idx = 1:bits_per_sample

            b = double(rand > 0.5);          % transmitted bit {0,1}

            if b == 1
                r = i1   + s1  * randn;      % signal + shot + thermal
            else
                r = 0    + sig0 * randn;      % thermal only
            end

            b_hat = double(r > tau);         % threshold decision

            if b_hat ~= b
                total_errors = total_errors + 1;
            end
            total_bits = total_bits + 1;
        end

        if mod(k, 1000) == 0
            fprintf('  k = %4d/%4d   running BER = %.3e\n', ...
                    k, N, total_errors/total_bits);
        end
    end

    BER_MC(p_idx) = total_errors / total_bits;

    fprintf('  MC  BER = %.4e  |  BER_new = %.4e  |  BER_old = %.4e\n\n', ...
            BER_MC(p_idx), BER_new(p_idx), BER_old(p_idx));
end

% -------------------------------------------------------------------------
% 5) SUMMARY TABLE
% -------------------------------------------------------------------------
fprintf('\n%s\n', repmat('=', 1, 72));
fprintf('%-10s  %-10s  %-14s  %-14s  %-14s\n', ...
        'P0 (dBm)', 'SNR (dB)', 'BER_MC', 'BER_new', 'BER_old');
fprintf('%s\n', repmat('-', 1, 72));
for i = 1:num_P0
    fprintf('%-10.1f  %-10.1f  %-14.4e  %-14.4e  %-14.4e\n', ...
            P0_dBm(i), SNR_dB_axis(i), BER_MC(i), BER_new(i), BER_old(i));
end
fprintf('%s\n', repmat('=', 1, 72));

% -------------------------------------------------------------------------
% 6) PLOTS
% -------------------------------------------------------------------------

% --- Plot 1: BER vs SNR (dB) ---
figure('Name', 'BER vs SNR — physical noise model');
semilogy(SNR_dB_axis, BER_MC,  'ko-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
semilogy(SNR_dB_axis, BER_new, 'bs--','LineWidth', 1.5, 'MarkerSize', 6);
semilogy(SNR_dB_axis, BER_old, 'r^--','LineWidth', 1.5, 'MarkerSize', 6);
yline(1e-9, 'k:', 'BER = 10^{-9}', 'LabelOrientation','horizontal', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)  [derived from P_0 and detector physics]');
ylabel('BER');
title('BER vs SNR — OOK with pointing jitter + physical detector noise');
legend('Monte Carlo', ...
       'Analytical: averaged over fading (BER_{new})', ...
       'Analytical: mean-I only (BER_{old})', ...
       'Location', 'southwest');

% --- Plot 2: BER vs received power (dBm) — matches Alice's plot style ---
figure('Name', 'BER vs Received Power');
semilogy(P0_dBm, BER_MC,  'ko-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
semilogy(P0_dBm, BER_new, 'bs--','LineWidth', 1.5, 'MarkerSize', 6);
semilogy(P0_dBm, BER_old, 'r^--','LineWidth', 1.5, 'MarkerSize', 6);
yline(1e-9, 'r--', 'BER = 10^{-9}', 'LabelOrientation','horizontal', 'LineWidth', 1.2);
grid on;
xlabel('Average received power P_0 (dBm)');
ylabel('BER');
title('BER vs Received Power — OOK with pointing jitter + physical detector noise');
legend('Monte Carlo', ...
       'Analytical: averaged over fading (BER_{new})', ...
       'Analytical: mean-I only (BER_{old})', ...
       'Location', 'southwest');