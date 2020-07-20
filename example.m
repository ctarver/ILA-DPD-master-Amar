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
shifts = [-0.05 + 0.05j, 0.05, 0.05, -0.1 - 0.05j, 0.05, 0.05, -0.1 - 0.05j, 0.05, 0.05]'; % A somewhat odd way to move coefficients in a square

% Move coefficients in a square grid, from left to right going from row by row
% [1]      [2]       [3] 
% [4] [5 (original)] [6] 
% [7]      [8]       [9]

shifts = [-0.05 + 0.05j, 0.05, 0.05, -0.1 - 0.05j, 0.05, 0.05, -0.1 - 0.05j, 0.05, 0.05]'; 
% all_coeffs_cell holds all coefficients along with the 'after' values to
% be plotted later
all_coeffs_cell = {};
prec_iter = 0;
precision = 4;
% If the minimum value is at the previously found coefficient, decrease the grid by half.
% 'precision' is a limit as to how many times this grid is cut in half
precision = 2; 

while prec_iter < precision
    % all_coeffs holds just the current grid 
    all_coeffs = [];

    for k1 = 1:length(shifts)
        for b=1:length(dpd.coeffs)
            % Loop through dpd.coeffs and apply a shift
            dpd.coeffs(b, 1) = dpd.coeffs(b, 1) + shifts(k1, 1);
        end


        % Transmit
        [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
        before = w_out_dpd.measure_all_powers;
        after = w_dpd.measure_all_powers;


        % Make copies
        inter_coeffs = dpd.coeffs;
        after_coeffs = after;


        % Append the leftmost value of after_coeff to the vector containing
        % the coefficients
        inter_coeffs = [inter_coeffs; after_coeffs(1,1)];

        % Now append this vector containing the coefficients with their
        % corresponding after value to all_coeffs. (This vector is reset in
        % every loop and is only used to compare 'after' values)
        all_coeffs = [all_coeffs, inter_coeffs];

        % Also append this vector to a cell vector used to make the plot at
        % the end
        all_coeffs_cell = [all_coeffs_cell, num2cell(inter_coeffs)];
    end

    disp(all_coeffs)

    % Find the minimum after value in the grid of 9 coefficient/after
    % vectors in all_coeffs
    [m, I] = min(all_coeffs(3, :));
    small = all_coeffs(:, I);
    disp(I)

    % If the coefficients corresponding to the smallest 'after' value is
    % the 5th/center value (aka the 'current' coefficient), then divide all
    % the shifts by half, increase prec_iter, and rerun the loop
    if I == 5
        for d=1:length(shifts)
            shifts(d) = shifts(d)/2;
        end
        prec_iter = prec_iter + 1;
        fprintf ("prec iter is: %d\n", prec_iter)
    end

    % Change current coefficient to the ones corresponding to the smallest
    % 'after' value
    for c=1:length(dpd.coeffs)
        dpd.coeffs(c,1) = all_coeffs(c, I);
    end
end
disp(all_coeffs)

% Plotting code

x_axis = {};
y_axis = {};
z_axis = {};
for order_iter = 1:length(dpd.coeffs)
    thisfig = figure();
    x_axis = {};
    y_axis = {};
    z_axis = {};
    for shift_iter = 1:length(all_coeffs_cell)
        x_axis = [x_axis, real(all_coeffs_cell{order_iter, shift_iter})];
        y_axis = [y_axis, imag(all_coeffs_cell{order_iter, shift_iter})];
        z_axis = [z_axis, all_coeffs_cell{length(dpd.coeffs) + 1, shift_iter}];
    end
    x_axis = cell2mat(x_axis);
    y_axis = cell2mat(y_axis);
    z_axis = cell2mat(z_axis);
    disp(x_axis);
    disp(y_axis);
    disp(z_axis);
    plot3(x_axis, y_axis, z_axis, 'o');
    xlabel('Real')
    ylabel('Complex')
    zlabel('after')
    if order_iter == 1
        title('1st coefficient');
    elseif order_iter == 2
        title('2nd coefficient');
    elseif order_iter == 3
        title('3rd coefficient');
    else
        title(strcat(num2str(order_iter), 'th coefficient'));
    end
    grid on
end
%% Plot
        w_out_dpd.plot_psd;
        w_dpd.plot_psd;
        dpd.plot_history;
