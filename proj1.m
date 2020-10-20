% GENERAL INFO: 
    % ->Code takes ~ 30 seconds to compile and get executed. 
    % ->In command window, the custom input (myMessage) will be displayed
    %   and 15 decoded signals under different SNR values will be displayed
    % ->Letter error rate and Word error rate plots will be generated
    % ->In "Decoding custom input message" section,
    %   There is already a default custom input containing all the character
    %   set,you can edit the custom input by editing the string array 
    %   in line no. 78    
    % ->IMPORTANT: YOU MUST EVALUATE INITIALIZATION SECTION (just defined below) to
    %   evaluate any of the main code sections i.e. LETTER ERROR RATE, WORD 
    %   ERROR RATE and DECODING CUSTOM INPUT MESSAGE, otherwise errors
    %   occur (we need the initialized data and variables for the entire
    %   code)
    % -> It is better to run the entire code in one go.
%% Initialization
% morse code mapping(bit serialization)
% making variables with the names corresponding to morse code value

% alphabets
A = [1 0 1 1 1];
B = [1 1 1 0 1 0 1 0 1];
C = [1 1 1 0 1 0 1 1 1 0 1];
D = [1 1 1 0 1 0 1];
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
% blank space between words
space = [0 0 0 0 0 0 0];

% Entire character set 
charSet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'];

% mapping the character set to their morse code bit equivalent using hash map / dictionary
keySet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ' '};
valueSet = {A, B, C, D, E, F, G, H, I , J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, n1, n2, n3, n4, n5, n6, n7, n8, n9, n0, space};
% Creating hash map/dictionary
map = containers.Map(keySet,valueSet);

% spacing between letters for morse
letterSpace = 3; 
% spacing between words for morse
wordSpace = 7; 

%% Decoding custom input message

%*****TRANSMITTER PART**************************************************
%*****ENCODER
% Sending a custom input, and decoding it.
myMessage = 'the quick brown fox jumped over the lazy dog 456127893'
% splitting the msg and now myMessage is a cell array
myMessage = split(myMessage);
% ->changing to uppercase (our morse code mapping is done only for upper
%   case characters)
myMessage= upper(myMessage);
% N is size of myMessage
N = length(myMessage);

% ->empty cell array initialized to store the morse bit equivalent of the
%   myMessage stored in bits cell array
bits = {};

% ->Converting each word in myWords cell array to their equivalent bits and
%   storing them in myWords_in_bits cell array
for i = 1:N
    for j =1:length(myMessage{i})
        % getting each character's morse code bits from myMessage
        char = map(myMessage{i}(j));      
        % appending the morse code bits to bits cell array
        bits{end+1} = char;
        if(j ~= length(myMessage{i}))
            % appending 3 zeroes after getting each character
            % We don't do it for last character thus the if condition
            bits{end+1} = zeros(1, letterSpace);
        end
    end
    
    if(i ~= N)
        % appending 7 zeroes after getting each word
        % We don't do it for last word thus the if condition
        bits{end+1} = zeros(1, wordSpace);
    end
end

%*****BPSK MODULATION
% here, mapping 1 -> -1 and 0 -> 1
myMessage_mod = BPSK_mod(bits);


%*****RECEIVER PART**************************************************
% a cell array is initilized to store 15 different decoded msgs under
% different SNR varying from 1-15
decoded_signals = {};

% we get decoded message under different SNR values
for SNR =1:1:15
    % if decoded bits dont match anything from char set, we show invalid
    invalid_char = '?';
    
    %*****AWGN NOISE ADDED
    % awgn noise added to bpsk modulated signal to get received signal
    received_signal = cellfun(@(x) awgn(x,SNR),myMessage_mod,'un',0);
    
    %*****BPSK DEMODULATION
    % converting 1-> 0 and -1 -> 1
    received_signal_in_bits = cellfun(@(x) fix(x*0 + (x<0)), received_signal, 'un', 0);
    % cell array to store decoded message
    decoded_message = {};

    %*****DECODER
    % we iterate over the cell array containing our received bits in each
    % cell, and decode the message
    for i = 1:length(received_signal_in_bits)
        % hasMatch is a boolean variable that is true when received bits matches any of
        % the morse code bits belonging to our valueSet(a set where all
        % our morse code mappings exist)
        hasMatch = any(cellfun(@isequal, valueSet, repmat({received_signal_in_bits{i}}, size(valueSet))));
        if mod(i,2) == 1 % odd indexes of cell array only has character bits
            if hasMatch == 1 % if set of bits belong to any morse code bits
                for j =1:length(keySet)
                    if isequal(map(keySet{j}), received_signal_in_bits{i})
                        % if received bits match any of the valueSet, we
                        % append the letter to decoded msg
                        decoded_message{end+1} = keySet{j};
                        break;
                    end
                end
            else
                % if not matching, we put invalid char
                decoded_message{end+1} = invalid_char;
            end
        else % even indexes of cell array have only letter spaces and word spaces
            if isequal(received_signal_in_bits{i}, [0 0 0 0 0 0 0])
                % appending a space
                decoded_message{end+1} = ' ';
            else
                if ~isequal(received_signal_in_bits{i}, [0 0 0])
                    % appending an invalid char
                    decoded_message{end+1} = invalid_char;
                end
            end
        end 
    end
    final_decoded_message = cell2mat(decoded_message); % the final message will be in capitals
    % storing the decoded message array as a cell in the decoded_signals
    % cell array
    decoded_signals{end+1} = final_decoded_message;
end
% displaying the cell array 
celldisp(decoded_signals);
%% LETTER ERROR RATE 

%*****TRANSMITTER PART**************************************************
%*****ENCODER
% empty cell array initialized to store 100 letters, where each cell is a letter
letters = {}; 

% Creating a cell array of 100 cells where each cell is a random letter
for i = 1:100
    % selecting random element from charSet
    pos = randi(length(charSet));
    letter = charSet(pos);
    
    % appending random element to letters array
    letters{end+1} = letter;
end

% ->empty cell array initialized to store the morse bit equivalent of the
%    earlier random letters stored in letters cell array
letters_in_bits = {};

% ->Converting each letter in letters cell array to their equivalent bits and
%   storing them in letters_in_bits cell array
for i = 1:length(letters)
    % ->going through each cell in letters cell array, and mapping them to
    %   their morse code bit equivalent
    letter_in_bit = map(letters{i});
    
    % ->appending the mapped morse code bit equivalent to letters_in__bits
    %   cell array
    letters_in_bits{end+1} = letter_in_bit;
end

%*****BPSK MODULATION
% here, mapping 1 -> -1 and 0 -> 1
letters_mod = BPSK_mod(letters_in_bits);

%*****ERROR RATE
% Make Letter Error Rate Graph
% The receiever blocks are present in makeErrorRateGraph Function
LER = makeErrorRateGraph(letters_in_bits, letters_mod, 1, 'Letter Error Rate');


%% WORD ERROR RATE

% ENCODER

% myWords is cell array of 20 words, here each cell is a random word
myWords = {'ametropia', 'syssitia', 'fever', 'capitule', 'tabellions', 'tutiorism', 'cannon', 'endosteum', 'apepsy', 'tupu', 'machzorim', 'unscale', 'presidia', 'cnicin', 'clownism', 'chittah', 'zaddiks', 'cuvy', 'quarry', 'enones'};

% changing to uppercase (our morse code mapping is done only for upper
% case characters)
myWords = cellfun(@(x) upper(x), myWords, 'un', 0);

% ->empty cell array initialized to store the morse bit equivalent of the
%   earlier random words stored in myWords cell array
myWords_in_bits = {};

% ->Converting each word in myWords cell array to their equivalent bits and
%   storing them in myWords_in_bits cell array
for i = 1:length(myWords)
    % empty array to store letters to make a word
    word = [];
    for j =1:length(myWords{i})
        % getting each letter from a particular random word
        letter = map(myWords{i}(j));     
        % appending the letter to word array 
        word = [word letter];
        if(j ~= length(myWords{i})) 
            % padding with 3 zeroes, (space between letters of a word is 3)
            % Also we don't want zeros at the end (thus the if condition)
            word = [word zeros(1, letterSpace)];
        end
    end
    % each word is appended to myWords_in_bits cell array
    myWords_in_bits{end+1} = word;
end

% BPSK MODULATION
% here, mapping 1 -> -1 and 0 -> 1
myWords_mod = BPSK_mod(myWords_in_bits);

% Make Word Error Rate Graph
% The receiever blocks are present in makeErrorRateGraph Function
WER = makeErrorRateGraph(myWords_in_bits, myWords_mod, 2, 'Word Error Rate');


%% FUNCTION DEFINITIONS

% counting num of errors, errors occur when cells of A and B are inequal
function [noe] = numerr(A, B)
    noe = 0;
    for i = 1:length(A)
        if ~isequal(A{i}, B{i})
            noe = noe + 1;
        end
    end
end

% makes error rate graphs i.e. error rate vs SNR
function [final_error_rate_list] = makeErrorRateGraph(symbol_in_bits, signal_mod, nth_graph, y_label)
    %*****RECEIVER PART**************************************************
    % this is used for plotting y axis
    final_error_rate_list = [];
    % size of modulated signal
    N = length(signal_mod);  
    
    % we take average of 1000 iterations for each SNR value, and get the
    % error rate for each SNR in final_error_rate_list array 
    for  SNR = 1:1:20
        error_rate_list = [];
        for i=1:1000
            %*****AWGN NOISE ADDED
            % adding awgn to signal sent via transmitter
            received_signal = cellfun(@(x) awgn(x,SNR),signal_mod,'un',0);
            
            %*****BPSK DEMODULATION
            % mapping -1 -> 1 and 1 -> 0
            received_signal_in_bits = cellfun(@(x) fix(x*0 + (x<0)), received_signal, 'un', 0);
            
            %*****ERROR RATE
            % getting num of errors between received signal and transmitted
            % signal
            noe = numerr(symbol_in_bits, received_signal_in_bits);
            % calculating error rate
            error_rate= noe/N;
            % appending error rate into error_rate_list array
            error_rate_list= [error_rate_list error_rate];  
        end
        % appending the average of 1000 iterations of each snr value to final_error_rate_list 
        final_error_rate_list = [final_error_rate_list mean(error_rate_list)];
    end
    
    % plotting
    SNR = 1:1:20;
    subplot(3,1,nth_graph);
    plot(SNR, final_error_rate_list, 'r');
    title(strcat(y_label,' vs SNR'));
    ylabel(y_label);
    xlabel('SNR');
end

% bpsk modulation is done to signal
function [mod] = BPSK_mod(signal)
    mod= cellfun(@(x) x*(-2),signal,'un',0);
    mod= cellfun(@(x) x+(1),mod,'un',0);
end