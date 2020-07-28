classdef PA_NN_Model < handle
    %PAModel Uses a NN to create a forward model of a PA.
    
    properties
        Class = 'PA_NN_Model'
    end
    
    properties
        n_epochs
        n_neurons
        n_hidden_layers
        memory_depth
        activation_function
        loss_function
        optimizer
        rnn % enable rnn
        delay % rnn delay
        use_dc_term % use a dc term
        learning_rate % How much of the new iteration to use vs previous iteration. Should be in (0, 1]
        net
        tr % for performance plot
        physical_atttenuation
        load_nn % to load existent net or not
        model % model name to be loaded
        history % Holds onto the coeffs used at each iteration
        sampling_rate
    end
    
    methods
        function obj = PA_NN_Model()
            obj.n_neurons = 12;
            obj.n_hidden_layers = 2;
            obj.memory_depth = 2;
            obj.activation_function = 'poslin';
            obj.loss_function = 'mse';
            obj.optimizer = 'trainlm';
            obj.rnn = 0;
            obj.delay = 1;
            obj.learning_rate = 0.01;
            obj.n_epochs = 1000;
            obj.load_nn = 0;
            obj.model = "";
            
            obj.use_dc_term = 1;
            obj.sampling_rate = 200e6;
            
            obj.run_setup();
        end
        
        function run_setup(obj)
            % setup neural network size
            nn_size = [];
            for i = 1:obj.n_hidden_layers
                nn_size = [nn_size, obj.n_neurons];
            end
            
            %net.numInputs = 2*obj.memory_depth;
            
            % initialize & specify optimizer
            if obj.rnn
                obj.net = layrecnet(obj.delay, nn_size, obj.optimizer);
                %             elseif obj.memory_depth > 1
                %                 net = timedelaynet(1:obj.memory_depth, nn_size, obj.optimizer);
            else
                obj.net = feedforwardnet(nn_size, obj.optimizer);
            end
            obj.net = init(obj.net);
            
            % set activation
            for i = 1:obj.n_hidden_layers
                obj.net.layers{i}.transferFcn = obj.activation_function;
            end
            
            % set loss function
            obj.net.performFcn = obj.loss_function;
            
            % set maximum epochs & learning rate
            obj.net.trainParam.epochs = obj.n_epochs;
            obj.net.trainParam.lr = obj.learning_rate;
        end
        
        function learn_model(obj, tx_signal, rx_signal)
            % learn PA NN Model considering memory effect
            % TODO: Try to delay manually
            tx_train = [real(tx_signal.data) imag(tx_signal.data)]';
            rx_train = [real(rx_signal.data) imag(rx_signal.data)]';
            tx_tmp = tx_train';
            for depth = 1:obj.memory_depth-1
                delayed_real = real(fft(delayseq(ifft(tx_signal.data), depth)));
                delayed_imag = imag(fft(delayseq(ifft(tx_signal.data), depth)));
                tx_tmp = [tx_tmp delayed_real delayed_imag];
            end
            tx_tmp = tx_tmp';
            
            obj.net = train(obj.net, tx_tmp, rx_train, 'useParallel','yes','useGPU','yes','showResources','yes');
        end
        
        function [out, raw] = transmit(obj, in, processing_flag)
            pa_in = in;
            if strcmp(obj.Class, 'webRF')
                pa_in = obj.pad_zeros(pa_in, 200);
                pa_in = obj.make_multiple_copies(pa_in, 4);
            end
            out = 0; % dummy. 
            raw = use_pa(obj, pa_in); % VSG uses the name to save to device
        end
        
        function out = use_pa(obj, input_signal)
            % main method to use the learned PA Model.
            input = [real(input_signal) imag(input_signal)]';
            input_tmp = input';
            for depth = 1:obj.memory_depth-1
                delayed_real = real(fft(delayseq(ifft(input_signal), depth)));
                delayed_imag = imag(fft(delayseq(ifft(input_signal), depth)));
                input_tmp = [input_tmp delayed_real delayed_imag];
            end
            input_tmp = input_tmp';
            
            output = obj.net(input_tmp, 'useParallel','yes','useGPU','yes','showResources','yes');
            out = 1j * output(2,:)' + output(1,:)';
        end
        
        function save(obj)
            % save neural net
            path = sprintf("%d_neurons_%d_layers_%d_memorydepth_%s_%s_nn", obj.n_neurons, ...
                obj.n_hidden_layers, obj.memory_depth, obj.loss_function, ...
                obj.activation_function);
            my_net = obj.net;
            save(path, "my_net");
        end
        
        function load(obj)
            % load neural net
            obj.net = load(obj.model,'net').net;
        end
        
        function out = rx_only(obj)
            out = obj.make_thermal_noise(10000);
        end
        
        function out = make_thermal_noise(obj, n_samples)
            noise_power = -174 + 10*log10(obj.sampling_rate);
            noise_vector = randn(n_samples, 1) + 1i*rand(n_samples, 1);
            noise_signal = Signal(noise_vector, obj.sampling_rate);
            noise_signal.normalize_to_this_rms(noise_power);
            out = noise_signal;
        end
        
        function plot_history(obj)
            % plot_history. Plots how the magnitude of the DPD coeffs
            % evolved over each iteration.
            figure(55);
            hold on;
            plot(abs(obj.history'));
            title('History for DPD Coeffs Learning');
            xlabel('Iteration Number');
            ylabel('abs(coeffs)');
            grid on;
        end
    end
end
