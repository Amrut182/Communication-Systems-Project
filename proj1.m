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
space = [0 0 0 0 0 0 0];

% Entire character set 
charSet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'];

% making hash map
keySet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ' '};
valueSet = {A, B, C, D, E, F, G, H, I , J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, n1, n2, n3, n4, n5, n6, n7, n8, n9, n0, space};
map = containers.Map(keySet,valueSet);
letterSpace = 3; % spacing between letters
wordSpace = 7; % spacing between words

% LETTER ERROR RATE
% cell array where each cell is a letter
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
    letter_in_bit = map(letters{i});
    letters_in_bits{end+1} = letter_in_bit;
end

% bpsk modulation
letters_mod = cellfun(@(x) x*(-2),letters_in_bits,'un',0);
letters_mod = cellfun(@(x) x+(1),letters_mod,'un',0);

% Make Letter Error Rate Graph
LER = makeErrorRateGraph(letters_in_bits, letters_mod, 1, 'Letter Error Rate');

%------------------
% WORD ERROR RATE
% cell array where each cell is a word
myWords = {'ametropia', 'syssitia', 'fever', 'capitule', 'tabellions', 'tutiorism', 'cannon', 'endosteum', 'apepsy', 'tupu', 'machzorim', 'unscale', 'presidia', 'cnicin', 'clownism', 'chittah', 'zaddiks', 'cuvy', 'quarry', 'enones'};
N_words = length(myWords)

% Always change to uppercase
myWords = cellfun(@(x) upper(x), myWords, 'un', 0)

% mapping words to bits
myWords_in_bits = {};
for i = 1:length(myWords)
    word = [];
    for j =1:length(myWords{i})
        letter = map(myWords{i}(j));      
        word = [word letter];
        if(j ~= length(myWords{i}))
            word = [word zeros(1, letterSpace)];
        end
    end
    myWords_in_bits{end+1} = word;
end

% bpsk modulation
myWords_mod= cellfun(@(x) x*(-2),myWords_in_bits,'un',0);
myWords_mod= cellfun(@(x) x+(1),myWords_mod,'un',0);

% Make Word Error Rate Graph
WER = makeErrorRateGraph(myWords_in_bits, myWords_mod, 2, 'Word Error Rate');

%------------
% Sending a custom input, and decoding it.
myMessage = 'the quick brown fox jumped over the lazy dog'
myMessage = split(myMessage);
myMessage= upper(myMessage);
N = length(myMessage);

bits = {};
for i = 1:N
    for j =1:length(myMessage{i})
        char = map(myMessage{i}(j));      
        bits{end+1} = char;
        if(j ~= length(myMessage{i}))
            bits{end+1} = zeros(1, letterSpace);
        end
    end
    if(i ~= N)
        bits{end+1} = zeros(1, wordSpace);
    end
end

% bpsk modulation
myMessage_mod= cellfun(@(x) x*(-2),bits,'un',0);
myMessage_mod= cellfun(@(x) x+(1),myMessage_mod,'un',0);

% NOISE
SNR = 5;
invalid_char = '?';
received_signal = cellfun(@(x) awgn(x,SNR),myMessage_mod,'un',0);
decoded_signal_in_bits = cellfun(@(x) fix(x*0 + (x<0)), received_signal, 'un', 0);
decoded_signal = [];
final_msg = [];
temp = {};

for i = 1:length(decoded_signal_in_bits)
    hasMatch = any(cellfun(@isequal, valueSet, repmat({decoded_signal_in_bits{i}}, size(valueSet))));
    if mod(i,2) == 1 % odd indexes have characters
        if hasMatch == 1
            for j =1:length(keySet)
                if isequal(map(keySet{j}), decoded_signal_in_bits{i})
                    temp{end+1} = keySet{j};
                    break;
                end
            end
        else
            temp{end+1} = invalid_char;
        end
    else % even indexes have only letter spaces and word spaces
        if isequal(decoded_signal_in_bits{i}, [0 0 0 0 0 0 0])
            temp{end+1} = ' ';
        else
            if ~isequal(decoded_signal_in_bits{i}, [0 0 0])
                temp{end+1} = invalid_char;
            end
        end
    end 
end

final_decoded_message = cell2mat(temp) % the final message will be in capitals

% FUNCTION DEFINITIONS
function [noe] = numerr(A, B)
    noe = 0;
    for i = 1:length(A)
        if ~isequal(A{i}, B{i})
            noe = noe + 1;
        end
    end
end

function [final_error_rate_list] = makeErrorRateGraph(symbol_in_bits, signal_mod, nth_graph, y_label)
    final_error_rate_list = [];
    N = length(signal_mod);  
    
    for  SNR = 1:1:20
        error_rate_list = [];
        for i=1:1000
            received_signal = cellfun(@(x) awgn(x,SNR),signal_mod,'un',0);
            decoded_signal = cellfun(@(x) fix(x*0 + (x<0)), received_signal, 'un', 0);
            noe = numerr(symbol_in_bits, decoded_signal);
            error_rate= noe/N;
            error_rate_list= [error_rate_list error_rate];  
        end
        final_error_rate_list = [final_error_rate_list mean(error_rate_list)];
    end
    
    SNR = 1:1:20;
    subplot(3,1,nth_graph);
    plot(SNR, final_error_rate_list, 'r');
    title(strcat(y_label,' vs SNR'));
    ylabel(y_label);
    xlabel('SNR');
end