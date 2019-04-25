clear;
clc;
format long;
%This code takes 4 sybols,alpahabets with uniform probabilities pr.
%N is number of symbols to be encoded, S is random genrated symbols.
alphabet=[0,1,2,3,4,5,6,7,8,9,10];
M = length(alphabet);
pr = ones(1,M);
C = [0,cumsum(pr)];
N = 100000;
S = round(rand(1,N)*max(alphabet));
ac_encoded_stream = zeros(0, 0, 'uint8');
ac_encoded_stream = arithmetic_encoder(N, S, M, C, ac_encoded_stream);
pr = ones(1,M);
C = [0,cumsum(pr)];
ac_encoded_stream = double(ac_encoded_stream);
symbol_stream = [];
symbol_stream = arithmetic_decoder(N, M, C,ac_encoded_stream, symbol_stream);
num_bits = 4*length(S); %4 bits required to represent 11 symbols.
fprintf(' \n total bits without compression =  %d',num_bits);
fprintf(' \n total bits after compression =  %d',length(ac_encoded_stream));
compression = (1 - (length(ac_encoded_stream)/num_bits)) * 100;
fprintf('\n compression percentage =  %.2f',compression);
if((S-symbol_stream) == 0)
fprintf(' \n no error in encoding and decoding ');
end

function d = arithmetic_encoder(N, S, M, C, d)
b = 0;
l = 1;
t = 0;
for k=1:N
    s_k = S(k);
    [b, l] = interval_update(s_k, b, l, M, C);
    if (b >= 1)
        b = b - 1;
        d = propagate_carry(t,d);
    end
    if (l <= 0.5)
        [b, l, t, d] = encoder_renormalization(b, l, t, d);
    end
    [C]  = update_distribution(s_k,M,C);
end
[d] = code_value_selection(b, t, d);
end
function [b,l] = interval_update(s, b, l, M, C)
g = l/C(M+1);
if (s == M-1)
    y = b+l;
else
    y = b + g*C(s+ 1 +1);
end
b = b+g*C(s +1);
l = y - b;
end
function [C]  = update_distribution(s,M,C)
 for m = s+1 : M
   C(m+1) = C(m+1) + 1;
 end
end
function [d] = propagate_carry(t, d)
n = t;
while(d(n) == 1)
    d(n)= 0;
    n = n - 1;
end
d(n) = 1;
end
function [b, l, t, d] = encoder_renormalization(b, l, t, d)
while (l <= 0.5)
    t = t+1;
    l = 2*l;
    if (b >= 0.5)
        d(t) = 1;
        b = 2*(b-0.5);
    else
        d(t) = 0;
        b = 2*b;
    end
end
end
function [d] = code_value_selection(b, t, d)
t = t + 1;
if (b <= 0.5)
    d(t) = 1;
else
    d(t) = 0;
    d = propagate_carry(t-1, d);
end
end
function S = arithmetic_decoder(N, M, C, d, S)
b = 0;
l = 1;
% assume P < # bits of in significant/matissa
%P = length(d);
%P = min(length(d), 52);
P = 50;

v = 0;
for n=1:P
    if (n <= length(d))
        v = v + d(n)/(2^n);
    end
end
t = P;
for k=1:N
    [s_k, l, b] = interval_selection(v, b, l, M, C);
    S(k) = s_k;
    if (b>= 1.0) 
        b = b - 1;
        v = v - 1;
    end
    if (l <= 0.5)
        [v, b, l, t] = decoder_renormalization(v, b, l, t, d, P);
    end
    [C]  = update_distribution(s_k,M,C);
end
end
function [s, l, b] = interval_selection(v, b, l, M, C)
s = M-1;
g = l/C(M+1);
x = b + g * C(M);
y = b + l;
while (x > v)
    s = s - 1;
    y = x;
    x = b + g * C(s+1);
end
b = x;
l = y - b;
end
function [v, b, l, t] = decoder_renormalization(v, b, l, t, d, P)
while (l <= 0.5)
    if (b >= 0.5)
        b = 2*(b - 0.5);
        v = 2*(v - 0.5);
    else
        b = 2*b;
        v = 2*v;
    end
    t = t + 1;
    if (t <= length(d))
        v = v + d(t)/(2^P);
    end
    l = 2*l;
end
end