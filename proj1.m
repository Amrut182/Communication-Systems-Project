% morse code mapping(bit serialization)
% making variables with the names corresponding to morse code value
% alphabets
A = [1 0 1 1 1];
B = [1 1 1 0 1 0 1 0 1];
C = [1 1 1 0 1 0 1 1 1 0 1];
D = [1 1 1 0 1 1 1 0 1];
E = [1];
F = [1 0 1 0 1 1 1 0 1];
G = [1 1 1 0 1 1 1 0 1];
H = [1 0 1 0 1 0 1];
I = [1 0 1];
J = [1 0 1 1 1 0 1 1 1 0 1 1 1];
K = [1 1 1 0 1 0 1 1 1];
L = [1 0 1 1 1 0 1 0 1];
M = [1 1 1 0 1 1 1];
N = [1 1 1 0 1];
O = [1 1 1 0 1 1 1 0 1 1 1];
P = [1 0 1 1 1 0];
Q = [1 1 1 0 1 1 1 0 1 0 1 1 1];
R = [1 0 1 1 1 0 1];
S = [1 0 1 0 1];
T = [1 1 1];
U = [1 0 1 0 1 1 1];
V = [1 0 1 0 1 0 1 1 1];
W = [1 0 1 1 1 0 1 1 1];
X = [1 1 1 0 1 0 1 0 1 1 1];
Y = [1 1 1 0 1 0 1 1 1 0 1 1 1];
Z = [1 1 1 0 1 1 1 0 1 0 1];

% numbers
n1 = [1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 1];
n2 = [1 0 1 0 1 1 1 0 1 1 1 0 1 1 1];
n3 = [1 0 1 0 1 0 1 1 1 0 1 1 1];
n4 = [1 0 1 0 1 0 1 0 1 1 1];
n5 = [1 0 1 0 1 0 1 0 1];
n6 = [1 1 1 0 1 0 1 0 1 0 1];
n7 = [1 1 1 0 1 1 1 0 1 0 1 0 1];
n8 = [1 1 1 0 1 1 1 0 1 1 1 0 1 0 1];
n9 = [1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1];
n0 = [1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 1];

% Entire character set 
charSet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'];

% making hash map
keySet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'};
valueSet = {A, B, C, D, E, F, G, H, I , J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, n1, n2, n3, n4, n5, n6, n7, n8, n9, n0};
Dict = containers.Map(keySet,valueSet);

letterSpace = 3; % spacing between letters
wordSpace = 7; % spacing between words

myMessage = "My name is Amrut";

% start uncomment
% Letters array
letters = {}; % a cell array to store 100 letters
for i = 1:100
    % selecting random element from charSet
    pos = randi(length(charSet));
    letter = charSet(pos);
    % appending random element to letters array
    letters{end+1} = letter;
end

% converting letters to bits
letters_in_bits = {};
for i = 1:length(letters)
    letter_in_bit = Dict(letters{i});
    letters_in_bits{end+1} = letter_in_bit;
%     % padding with 3 zeroes(no padding after last letter, hence the if condition)
%     if(i ~= length(letters))
%         letters_in_bits{end+1} = zeros(1, letterSpace);
%     end
end

% NRZ encoding
letters_nrz = cellfun(@(x) x*(-2),letters_in_bits,'un',0);
letters_nrz = cellfun(@(x) x+(1),letters_nrz,'un',0);
% celldisp(letters_nrz);

% Noise
ber_final_sim = {};
N = length(letters_nrz);

for  SNR = 1:1:15
    ber_sim = {};
    for i=1:10
%         celldisp(letters_in_bits);
        received_signal = cellfun(@(x) awgn(x,SNR),letters_nrz,'un',0);
%         celldisp(received_signal);
        decoded_signal = cellfun(@(x) fix(x*0 + (x<0)), received_signal, 'un', 0);
%         celldisp(decoded_signal);
        noe = numerr(letters_in_bits, decoded_signal);
        ber_sim1 = noe/N
%         ber_sim = [ber_sim ber_sim1];  
    end
%     n = cell2mat(ber_sim);
%     m = mean(ber_sim);
%     ber_final_sim = [ber_final_sim m];
end

% SNRdB = 1:1:15;
% semilogy(SNRdB, ber_final_sim, 'r');


% AWGN noise addition
% end of uncomment

% %start trial
% myLetter = [A B C D E];
% N=1000000;
% m=randi([0 1],1,N);
% myLetter = m;
% myLetter1 = m;
% % % bpsk modulation, changing 1 -> -1, 0 -> 1
% % bpskModulator = comm.BPSKModulator;
% % modData = bpskModulator(myLetter');
% % bpsk_myLetter= real(modData');
% 
% myLetter = myLetter*(-2);
% myLetter = myLetter + 1;
% 
% % creating and showing bpsk modulated signal from a series of bits
% % [bpsk_myLetter, t] = makeSignal(bpsk_myLetter);
% % x = bpsk_myLetter; % x is transmitted signal
% % subplot(3, 1, 1);
% % plot(t, x,'r');
% 
% % Noise
% ber_sim = [];
% for  SNRdB = 1:1:15
% SNR=10^(SNRdB/10); %convert to normal scale
% sigma=sqrt(1/(2*SNR));
% % myLetter A = 10111
% r = myLetter + sigma.*randn(1, length(myLetter));
% m_cap = (r<0);
% % for i = 1:length(r)
% % if r(i)<0
% %     m_cap = 1;
% % else
% %     m_cap = 0;
% % end
% % end
% noe = sum(myLetter1~=m_cap);
% ber_sim1 = noe/length(myLetter);
% ber_sim = [ber_sim ber_sim1];  
% end
% SNRdB = 1:1:15;
% semilogy(SNRdB, ber_sim, 'r');
% %end trial

% function for making digital signals from bits
function [signal, t] = makeSignal(bits)
    N = size(bits, 2);
    S = 100;
    i = 1;
    t = 0 : 1/S : N;
    for j=1:length(t)
        if t(j) <= i
            signal(j) = bits(i);
        else 
            signal(j) = bits(i);
            i = i + 1;
        end
    end
end

function [noe] = numerr(A, B)
    noe = 0;
    for i = 1:length(A)
%         x = A{i}
%         y = B{i}
%         z = noe
        if ~isequal(A{i}, B{i})
            noe = noe + 1;
        end
    end
end
