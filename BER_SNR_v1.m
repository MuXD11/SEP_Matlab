%% MONTE CARLO BER — Multi-SNR version
% Computes Monte Carlo BER and analytical BER for 4 SNR levels.

% -------------------------------------------------------------------------
% 0) SIMULATION PARAMETERS  (must match the control-loop script)
% -------------------------------------------------------------------------
% crossover frequency [Hz] and sample frequency [Hz] come from previous sim
Tstep = 1/fs;
Tsim  = 10;             % s
N     = Tsim/Tstep;     % number of fading samples

% -------------------------------------------------------------------------
% 1) COMMUNICATION PARAMETERS
% -------------------------------------------------------------------------
bit_rate        = 1e8;               % 100 Mbps
bits_per_sample = round(bit_rate/fs);% bits per fading block  (= 20 000)
threshold       = 0.5;               % OOK mod decision threshold

% -------------------------------------------------------------------------
% 2) SNR LEVELS  (dB → linear → noise std)
% -------------------------------------------------------------------------
SNR_dB_array = [20, 15, 10, 5];          % four SNR points
num_snr      = numel(SNR_dB_array);

sigma_array  = zeros(1, num_snr);
for i = 1:num_snr
    SNR_lin        = 10^(SNR_dB_array(i)/10);
    sigma_array(i) = 1/sqrt(2*SNR_lin);  % noise std for AWGN + OOK
end

% -------------------------------------------------------------------------
% 3) Q-FUNCTION. This function indicates the statiscial probability that
% x>Z for a gaussian distribution
% -------------------------------------------------------------------------
Q = @(x) 0.5 * erfc(x / sqrt(2));

% -------------------------------------------------------------------------
% 4) TRIM I_CLOSED TO THE BER SIMULATION LENGTH
% -------------------------------------------------------------------------
I_sim = I_closed(1:N);   % N samples at fs = 5000 Hz  →  10 s

% -------------------------------------------------------------------------
% 5) PREALLOCATE RESULT ARRAYS
% -------------------------------------------------------------------------
BER_MC      = zeros(1, num_snr);
BER_old     = zeros(1, num_snr);   % classical formula: uses only mean(I)
BER_new     = zeros(1, num_snr);   % improved formula:  averages over fading

% -------------------------------------------------------------------------
% 6) MAIN LOOP OVER SNR LEVELS
% -------------------------------------------------------------------------
for snr_idx = 1:num_snr

    sigma_n  = sigma_array(snr_idx);
    SNR_lin  = 10^(SNR_dB_array(snr_idx)/10);

    % --- 6a) Analytical BER estimates --------------------------------

    % Classical: Ignores fading distribution, uses only the mean):
    BER_old(snr_idx) = 0.5 * erfc(sqrt(SNR_lin * mean(I_sim) / 2));

    % Improved Averages BER for every fading sample using:
    %   Pb(I_k) = 0.5*Q(threshold/sigma_n)           [bit=0 term]
    %           + 0.5*Q((sqrt(I_k)-threshold)/sigma_n)[bit=1 term]
    Pb_blocks          = 0.5*(Q(threshold/sigma_n)+ Q((sqrt(I_sim) - threshold)/sigma_n));
    BER_new(snr_idx)   = mean(Pb_blocks);

    % --- 6b) Monte Carlo BER -----------------------------------------
    total_errors = 0;
    total_bits   = 0;

    fprintf('\n[SNR = %2d dB | sigma_n = %.4f]\n', ...
            SNR_dB_array(snr_idx), sigma_n);

     % Iterate over blocks (each sample of I_k contains 20000 bits
    for k = 1:N

        I_k = I_sim(k);     % fading coefficient for this block

        for bit_idx = 1:bits_per_sample

            b       = double(rand > 0.5);       % transmitted bit  {0,1}
            s_faded = b * sqrt(I_k);            % OOK + amplitude fading
            r       = s_faded + sigma_n*randn;  % received sample

            b_hat   = double(r > threshold);    % hard decision

            if b_hat ~= b
                total_errors = total_errors + 1;
            end
            total_bits = total_bits + 1;
        end

        % Progress report every 10000 fading blocks
        if mod(k, 10000) == 0
            fprintf('  k = %4d / %4d   running BER = %.3e\n', ...
                    k, N, total_errors/total_bits);
        end
    end

    BER_MC(snr_idx) = total_errors / total_bits;

    fprintf('  --> MC BER  = %.4e\n', BER_MC(snr_idx));
    fprintf('  --> BER_old = %.4e\n', BER_old(snr_idx));
    fprintf('  --> BER_new = %.4e\n', BER_new(snr_idx));
    fprintf('  Total bits  = %.2e  |  Total errors = %d\n', ...
            total_bits, total_errors);
end



% -------------------------------------------------------------------------
% 7) SUMMARY TABLE STATISTICS
% -------------------------------------------------------------------------
fprintf('\n%s\n', repmat('=',1,65));
fprintf('%-10s  %-14s  %-14s  %-14s\n', ...
        'SNR (dB)', 'BER_MC', 'BER_old', 'BER_new');
fprintf('%s\n', repmat('-',1,65));
for i = 1:num_snr
    fprintf('%-10d  %-14.4e  %-14.4e  %-14.4e\n', ...
            SNR_dB_array(i), BER_MC(i), BER_old(i), BER_new(i));
end
fprintf('%s\n', repmat('=',1,65));

%% MONTE CARLO BER — Multi-SNR version
% Computes Monte Carlo BER and analytical BER for 4 SNR levels.

% -------------------------------------------------------------------------
% 0) SIMULATION PARAMETERS  (must match the control-loop script)
% -------------------------------------------------------------------------
Tstep = 1/fs;
Tsim  = 10;             % s
N     = Tsim/Tstep;     % number of fading samples

% -------------------------------------------------------------------------
% 1) COMMUNICATION PARAMETERS
% -------------------------------------------------------------------------
bit_rate        = 1e8;                % 100 Mbps
bits_per_sample = round(bit_rate/fs);% bits per fading block
threshold       = 0.5;               % OOK decision threshold

% -------------------------------------------------------------------------
% 2) SNR LEVELS  (FIXED ORDER: LOW → HIGH)
% -------------------------------------------------------------------------
SNR_dB_array = [5, 10, 15, 20];      % ascending order
num_snr      = numel(SNR_dB_array);

sigma_array  = zeros(1, num_snr);
for i = 1:num_snr
    SNR_lin        = 10^(SNR_dB_array(i)/10);
    sigma_array(i) = 1/sqrt(2*SNR_lin);
end

% -------------------------------------------------------------------------
% 3) Q-FUNCTION
% -------------------------------------------------------------------------
Q = @(x) 0.5 * erfc(x / sqrt(2));

% -------------------------------------------------------------------------
% 4) TRIM I_CLOSED
% -------------------------------------------------------------------------
I_sim = I_closed(1:N);

% -------------------------------------------------------------------------
% 5) PREALLOCATE
% -------------------------------------------------------------------------
BER_MC  = zeros(1, num_snr);
BER_old = zeros(1, num_snr);
BER_new = zeros(1, num_snr);

% -------------------------------------------------------------------------
% 6) MAIN LOOP OVER SNR
% -------------------------------------------------------------------------
for snr_idx = 1:num_snr

    sigma_n = sigma_array(snr_idx);
    SNR_lin = 10^(SNR_dB_array(snr_idx)/10);

    % --- Analytical BERs ---
    BER_old(snr_idx) = 0.5 * erfc(sqrt(SNR_lin * mean(I_sim) / 2));

    Pb_blocks = 0.5*(Q(threshold/sigma_n) + ...
                    Q((sqrt(I_sim) - threshold)/sigma_n));
    BER_new(snr_idx) = mean(Pb_blocks);

    % --- Monte Carlo ---
    total_errors = 0;
    total_bits   = 0;

    fprintf('\n[SNR = %2d dB | sigma_n = %.4f]\n', ...
            SNR_dB_array(snr_idx), sigma_n);

    for k = 1:N

        I_k = I_sim(k);

        for bit_idx = 1:bits_per_sample

            b       = double(rand > 0.5);
            s_faded = b * sqrt(I_k);
            r       = s_faded + sigma_n*randn;

            b_hat = double(r > threshold);

            if b_hat ~= b
                total_errors = total_errors + 1;
            end
            total_bits = total_bits + 1;
        end

        % Progress report
        if mod(k, 10000) == 0
            fprintf('  k = %4d / %4d   running BER = %.3e\n', ...
                    k, N, total_errors/total_bits);
        end
    end

    BER_MC(snr_idx) = total_errors / total_bits;

    fprintf('  --> MC BER  = %.4e\n', BER_MC(snr_idx));
    fprintf('  --> BER_old = %.4e\n', BER_old(snr_idx));
    fprintf('  --> BER_new = %.4e\n', BER_new(snr_idx));
    fprintf('  Total bits  = %.2e  |  Total errors = %d\n', ...
            total_bits, total_errors);
end

% -------------------------------------------------------------------------
% 7) SUMMARY TABLE
% -------------------------------------------------------------------------
fprintf('\n%s\n', repmat('=',1,65));
fprintf('%-10s  %-14s  %-14s  %-14s\n', ...
        'SNR (dB)', 'BER_MC', 'BER_old', 'BER_new');
fprintf('%s\n', repmat('-',1,65));
for i = 1:num_snr
    fprintf('%-10d  %-14.4e  %-14.4e  %-14.4e\n', ...
            SNR_dB_array(i), BER_MC(i), BER_old(i), BER_new(i));
end
fprintf('%s\n', repmat('=',1,65));

% -------------------------------------------------------------------------
% 8) BER vs SNR PLOT (FIXED)
% -------------------------------------------------------------------------
figure('Name','BER vs SNR');

semilogy(SNR_dB_array, BER_MC,  'ko-', 'LineWidth',1.8, 'MarkerSize',8); hold on;
semilogy(SNR_dB_array, BER_old, 'r^--','LineWidth',1.5, 'MarkerSize',7);
semilogy(SNR_dB_array, BER_new, 'bs--','LineWidth',1.5, 'MarkerSize',7);

grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR — OOK with pointing jitter (closed-loop)');

legend('Monte Carlo', ...
       'Analytical: mean-I only (BER_{old})', ...
       'Analytical: averaged over fading (BER_{new})', ...
       'Location','northeast');