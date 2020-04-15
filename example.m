clear; clc; close all;
%% Setup Everything

% Add the submodules to path
addpath(genpath('OFDM-Matlab'))
addpath(genpath('WARPLab-Matlab-Wrapper'))
addpath(genpath('Power-Amplifier-Model'))

rms_input = 0.50;

% Setup the PA simulator or TX board
PA_board = 'webRF'; % either 'WARP', 'webRF', or 'none'

switch PA_board
    case 'WARP'
        warp_params.nBoards = 1;         % Number of boards
        warp_params.RF_port  = 'A2B';    % Broadcast from RF A to RF B. Can also do 'B2A'
        board = WARP(warp_params);
        Fs = 40e6;    % WARP board sampling rate.
    case 'none'
        board = PowerAmplifier(7, 4);
        Fs = 40e6;    % WARP board sampling rate.
    case 'webRF'
        dbm_power = -22;
        board = webRF(dbm_power);
end

% Setup OFDM
ofdm_params.nSubcarriers = 1200;
ofdm_params.subcarrier_spacing = 15e3; % 15kHz subcarrier spacing
ofdm_params.constellation = 'QPSK';
ofdm_params.cp_length = 144; % Number of samples in cyclic prefix.
ofdm_params.nSymbols = 14;
modulator = OFDM(ofdm_params);

% Create TX Data
[tx_data, ~] = modulator.use;
tx_data = Signal(tx_data, modulator.sampling_rate, rms_input);
tx_data.upsample(board.sample_rate)

% Setup DPD
dpd_params.order = 3;
dpd_params.memory_depth = 1;
dpd_params.lag_depth = 0;  % 0 is a standard MP. >0 is GMP.
dpd_params.nIterations = 1;
dpd_params.learning_rate = 0.75;
dpd_params.learning_method = 'newton'; % Or 'ema' for exponential moving average.
dpd_params.use_even = false; 
dpd_params.use_conj = 0;    % Conjugate branch. Currently only set up for MP (lag = 0)
dpd_params.use_dc_term = 0; % Adds an additional term for DC
dpd = ILA_DPD(dpd_params);

%% Run Experiment
[~, w_out_dpd] = board.transmit(tx_data.data);
dpd.perform_learning(tx_data.data, board);
% Loop over coefficients here
shifts = {0.05, -0.1, 0.1j, -0.1j}
for k1 = 1:length(shifts)
    for b=1:length(dpd.coeffs)
        dpd.coeffs(b,1) = dpd.coeffs(b, 1) + shifts{k1};
        [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));

        before = w_out_dpd.measure_all_powers;
        after = w_dpd.measure_all_powers;

        %% Plot
        w_out_dpd.plot_psd;
        w_dpd.plot_psd;
        dpd.plot_history;
    end
end
%for b=1:length(dpd.coeffs)
%    dpd.coeffs(b,1) = dpd.coeffs(b,1) - 0.05j
%end


