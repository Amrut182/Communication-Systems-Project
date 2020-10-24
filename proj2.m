% GENERAL INFO: -> The program will take about a minute to run
%               -> 5 plots will be generated for 5 values of a, each plot
%               will have 3 simulations (represented by different colors in
%               graph) for bit error rate vs SNR
%               -> in command window, you can see transmitted signal,
%               modulated signal, received signal, demodulated signal and
%               decoded signal (for a = 5 and SNR = 30)
%               -> 4 out of the 5 plots are always coming as expected, 
%               (expected here means the 3 colors(3 simulations) are nearly coinciding), 
%               sometimes 1-2 plots are showing error for some reason which im
%               not able to figure out. worst case scenario, run program
%               again if none of the plots are coming like any of the plots in 
%               the image I have attached in the zip file (named plots.jpg)

clear,clc;

% *****TRANSMITTED SIGNAL
% Signal we transmit
transmittedSignal = randi([0 1],10 ,1)
N = length(transmittedSignal);

% taking a range of values for a (from 1 to 5)
for a = 1:1:5
    %*****MODULATION
    % mapping 0->[a, 0] and 1->[0, a]
    % the 1st column of modulated signal are x components, the 2nd column
    % of modulatedSignal are y components
    modulatedSignal = zeros(N,2);
    for i = 1:N
        if transmittedSignal(i) == 0
            modulatedSignal(i,1) = a;
            modulatedSignal(i,2) = 0;
        else
            modulatedSignal(i,1) = 0;
            modulatedSignal(i,2) = a;        
        end
    end

    % h = (h_real) + j*(h_imaginary), N(0, 1)
    h_real = randn(1, 1);
    h_imaginary = randn(1, 1);
    % h_real will multiply the corresponding x elements
    % h_imaginary will multiply the corresponding y elements
    % these are scalar quantities, as h is a scalar quantity
    
%     ploting for 3 simulations
    for j = 1:1:3
        final_error_rate_list = [];
        
        % we take average of 1000 iterations for each SNR value, and get the
        % error rate for each SNR in final_error_rate_list array 
        for  SNR = 1:0.2:30
            error_rate_list = [];
            for i=1:1000
                % *****RECEIVED SIGNAL 
                % adding noise and channel characteristics to our 
                % transmitted signal
                receivedSignal = zeros(N,2);
                % Multiply h_real with corresponding x components of
                % modulated signal
                receivedSignal(:,1) = h_real*modulatedSignal(:,1);
                % Multiply h_imaginary with corresponding y components of
                % modulated signal
                receivedSignal(:,2) = h_imaginary*modulatedSignal(:,2);

                % we know that for same snr value, awgn will give different
                % results
                % applying awgn to column containing x components
                receivedSignal(:,1) = awgn(receivedSignal(:, 1),SNR);
                % applying awgn to column containing y components
                receivedSignal(:,2) = awgn(receivedSignal(:, 2),SNR);
                % hence different noises of same SNR are added to x and y components by
                % doing the above 2 steps
                % here n = [n_real n_imaginary] so we are adding noise to x
                % and y components respectively via awgn 
                
                % *****DEMODULATOR
                % Demodulating the received signal
                demodulatedSignal = zeros(N, 2);
                for i = 1:length(receivedSignal(:,1))% loop through all rows
                % if energy of x component is lesser than energy of y
                % component, then it is demodulated as [0 a]
                % we know that the received signal is more likely to be [0
                % a] if the energy of 2nd component is greater than 1st
                % component
                 if (receivedSignal(i,1).^2) < (receivedSignal(i,2).^2) %checking energies
                    demodulatedSignal(i, 1) = 0;
                    demodulatedSignal(i, 2) = a;
                 % else demodulated as [a 0]
                 else
                    demodulatedSignal(i,1) = a;
                    demodulatedSignal(i, 2) = 0;
                 end
                end
                
                % *****DECODER
                % Decoding the demodulated signal
                decodedSignal = zeros(N, 1);
                for i = 1:length(demodulatedSignal(:,1))                    % loop through all rows
                    if (demodulatedSignal(i,1) == a)
                        decodedSignal(i) = 0;
                    else
                        decodedSignal(i) = 1;
                    end
                end
                
                % *****ERROR RATE
                % finding number of errors, (checking how many elements are
                % unequal)
                noe = sum(transmittedSignal ~= decodedSignal);
                % calculating error rate
                error_rate= noe/N;
                % appending error rate into error_rate_list array
                error_rate_list= [error_rate_list error_rate];  
            end
            % appending the average of 1000 iterations of each snr value to final_error_rate_list 
            final_error_rate_list = [final_error_rate_list mean(error_rate_list)];
        end

        % plotting
        SNR = 1:0.2:30;
        subplot(5,1,a);
        color_of_plots = ['r', 'b', 'g'];
        plot(SNR, final_error_rate_list, color_of_plots(j));
%         plot(SNR, final_error_rate_list, color_of_plots);
        hold on;
        title(['BER vs SNR (For a =', num2str(a),')']);
        ylabel('BER');
        xlabel('SNR');
    end
    hold off;
end

% Just to show an eg. of the transmittedSignal at intermediate stages, (to 
% be seen in command window), here a = 5, and SNR = 30 
(modulatedSignal)
(receivedSignal)
(demodulatedSignal)
(decodedSignal)