% Simpsons one third rule of integration


clc       % Clear command window
clear     % Clear variables
close all % Close figures

syms f(t)  % Define symbolic function
f(t) = 2000 * (log(140000 / (140000 - 2100 * t))) - 9.8 * t;

a = 8;   % Lower limit
b = 30;  % Upper limit
n = 300;   % Sub-intervals (even for Simpson's rule)

v1 = double(int(f, t, a, b)); % Exact integral
h = (b - a) / n;  % Step size

s1 = double(f(a)); % First term
s2 = double(f(b)); % Last term

s3 = 0; % Sum for even indices
for i = a + 2 * h : 2 * h : b - 2 * h
    s3 = s3 + f(i);
end
s3 = 2 * s3;

s4 = 0; % Sum for odd indices
for i = a + h : 2 * h : b
    s4 = s4 + f(i);
end
s4 = 4 * s4;

v2 = double((h / 3) * (s1 + s2 + s3 + s4)); % Simpson's 1/3 rule

% Display exact integral, numerical result, and error
val = sprintf("%.6f %.6f %.6f", v1, v2, v1 - v2);
disp(val);
