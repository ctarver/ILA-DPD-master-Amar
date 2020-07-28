clear; clc; close all;
%% Setup Everything

% Add the submodules to path
addpath(genpath('OFDM-Matlab'))
addpath(genpath('WARPLab-Matlab-Wrapper'))
addpath(genpath('Power-Amplifier-Model'))

rms_input = 0.50;

% Setup the PA simulator or TX board
PA_board = 'NN'; % either 'WARP', 'webRF', ''NN', or 'none'

switch PA_board
    case 'WARP'
        warp_params.nBoards = 1;         % Number of boards
        warp_params.RF_port  = 'A2B';    % Broadcast from RF A to RF B. Can also do 'B2A'
        board = WARP(warp_params);
        Fs = 40e6;    % WARP board sampling rate.
    case 'none'
        board = PowerAmplifier(7, 4);
        Fs = 40e6;    % WARP board sampling rate.
    case {'webRF', 'NN'}
        dbm_power = -24; % Originally -22
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

% If we are working with the NN, broadcast through the rfweblab to get
% training data, train the NN, then use that going forward instead of
% RFWebLeb
switch PA_board
    case 'NN'
        [~, web_rf_w_out_dpd] = board.transmit(tx_data.data);
        board = PA_NN_Model();
        board.learn_model(tx_data, web_rf_w_out_dpd)
end


%% Run Experiment
[~, w_out_dpd] = board.transmit(tx_data.data);
dpd.perform_learning(tx_data.data, board);

[~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
after = w_dpd.measure_all_powers;

after_values = [after(1,1)];
total_gradient = inf;

while total_gradient > 0.5 % Iterates until the sum of the gradient vector is less than 0.5
    % Shift in R^2 space
    del = 0.05 + 0.05i;
    coeff_length = length(dpd.coeffs);
    
    % Initialize an empty gradient_vector
    gradient_vector = zeros(coeff_length,1);
    
    for i = 1:coeff_length % Fills in each value in the gradient vector
        original_coeffs = dpd.coeffs;
        plus_delta = dpd.coeffs;
        minus_delta = dpd.coeffs;
        
        % Plus_delta is the shift in the positive direction in R^2 whereas
        % minus_delta is shift in the negative direction
        plus_delta(i, 1) = plus_delta(i, 1) + del;
        minus_delta(i, 1) = minus_delta(i, 1) - del;
        
        %Ideally I would just make a new dpd object, but since I'm not sure
        %about the specifics, here I'm just resetting the dpd.coeffs to be
        %whatever I want them to be, transmitting that, and then restoring
        %dpd.coeffs back to what they were originally. It's a bit rachet,
        %but I'm not entirely sure of a better way to do it
        dpd.coeffs = plus_delta;
        [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
        after = w_dpd.measure_all_powers;
        plus_delta = after(1,1);
        
        dpd.coeffs = minus_delta;
        [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
        after = w_dpd.measure_all_powers;
        minus_delta = after(1,1);
        
        dpd.coeffs = original_coeffs;
        gradient_vector(i, 1) = (plus_delta - minus_delta)/sqrt(.02);
    end
    
    % This creates the Hessian Matrix
    hess_dim = length(dpd.coeffs);
    hessian_grid = zeros(hess_dim, hess_dim);
    
    % Filling in the Hessian Matrix
    for i = 1:hess_dim
        for j = 1:hess_dim
            original_coeffs = dpd.coeffs;
            % Filling in the second derivatives along the diagonal
            if i == j
                plus_del = dpd.coeffs;
                zero_del = dpd.coeffs;
                minus_del = dpd.coeffs;
                plus_del(i, 1) = plus_del(i, 1) + del;
                minus_del(i, 1) = minus_del(i, 1) - del;
                
                % Once again, I don't know how to do this other than just
                % redefining dpd.coeffs, transmitting that, and then
                % resetting dpd.coeffs back to the original one
                dpd.coeffs = plus_del;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                plus_del_out = after(1,1);
                
                dpd.coeffs = zero_del;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                zero_del_out = after(1,1);
                
                dpd.coeffs = minus_del;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                minus_del_out = after(1,1);
                
                % Reset dpd.coeffs back to the original coeffs
                dpd.coeffs = original_coeffs;
                % Set the value in the Hessian Matrix
                hessian_grid(i, j) = (plus_del_out - 2 * zero_del_out + minus_del_out) / (del^2);
                
                % For all other values not on the diagonal
            else
                % If you look at the formula I wrote out, I have to pass
                % different values into the "function", which in this case
                % is just transmitting the dpd object and observing the
                % output spectrum power. This is a way I figured out how to
                % do it, but I'll see if it's possible to simplify it down
                % a bit.
                copy1 = dpd.coeffs;
                copy2 = dpd.coeffs;
                copy3 = dpd.coeffs;
                copy4 = dpd.coeffs;
                copy1(j, 1) = copy1(j, 1) + del;
                copy1(i, 1) = copy1(i, 1) + del;
                copy2(j, 1) = copy2(j, 1) - del;
                copy2(i, 1) = copy2(i, 1) + del;
                copy3(j, 1) = copy3(j, 1) + del;
                copy3(i, 1) = copy3(i, 1) - del;
                copy4(j, 1) = copy4(j, 1) - del;
                copy4(i, 1) = copy4(i, 1) - del;
                
                % Transmitting the four different vectors
                dpd.coeffs = copy1;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                copy1_out = after(1,1);
                
                dpd.coeffs = copy2;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                copy2_out = after(1,1);
                
                dpd.coeffs = copy3;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                copy3_out = after(1,1);
                
                dpd.coeffs = copy4;
                [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
                after = w_dpd.measure_all_powers;
                copy4_out = after(1,1);
                
                % Reset dpd.coeffs back to the original coeffs
                dpd.coeffs = original_coeffs;
                % Set the value in the Hessian Matrix
                hessian_grid(i, j) = (copy1_out - copy2_out - copy3_out - copy4_out) / (4 * (del ^ 2));
            end
        end
    end
    inv_hess = inv(hessian_grid);
    
    %     % See how much optimizer is changing the dpd.coeff vector by
    %     new_vect = dpd.coeffs - (inv_hess * gradient_vector);
    %     disp(new_vect - dpd.coeffs);
    
    %xk+1 = xk - [Hf(xk)]^-1 * gradient(f(xk))
    dpd.coeffs = dpd.coeffs - (inv_hess * gradient_vector);
    
    total_gradient = sum(gradient_vector);
    disp(total_gradient)
    
    % Add to the after_values vector to tell whether model is actually
    % optimizing
    [~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
    after_lvalue = w_dpd.measure_all_powers;
    after_values = [after_values, after_lvalue(1,1)];
    
end

disp(after_values);

%% Plot
w_out_dpd.plot_psd;
w_dpd.plot_psd;
dpd.plot_history;

