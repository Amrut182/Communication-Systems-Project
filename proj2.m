clear,clc;
% InitimodulatedSignallizing
transmittedSignal = randi([0 1],10 ,1)
N = length(transmittedSignal);

% ModulmodulatedSignaltion
a = 2;
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

k = 10;
sigma = k;
% Noise_vector = normrnd(0,sqrt(sigma),[N,1]);
for j = 1:1:3
    final_error_rate_list = [];
    
    for  SNR = 1:0.1:30
        error_rate_list = [];
        for i=1:1000
            receivedSignal = zeros(N,2);
            receivedSignal(:,1) = h_real*modulatedSignal(:,1);
            receivedSignal(:,2) = h_imaginary*modulatedSignal(:,2);
            receivedSignal = awgn(receivedSignal,SNR);

            Y = zeros(N, 1);
            for i = 1:length(receivedSignal(:,1))                    % loop through all rows
             if (receivedSignal(i,1).^2) >(receivedSignal(i,2).^2) %checking energies
                Y(i) = 0;
             else
                Y(i) = 1;
             end
            end

            noe = sum(transmittedSignal ~= Y);
            % calculating error rate
            error_rate= noe/N;
            % appending error rate into error_rate_list array
            error_rate_list= [error_rate_list error_rate];  
        end
        % appending the average of 1000 iterations of each snr value to final_error_rate_list 
        final_error_rate_list = [final_error_rate_list mean(error_rate_list)];
    end

    % plotting
    SNR = 1:0.1:30;
    subplot(3,1,j);
    plot(SNR, final_error_rate_list, 'r');
    title(('BER vs SNR'));
    ylabel('BER');
    xlabel('SNR');
end


 

