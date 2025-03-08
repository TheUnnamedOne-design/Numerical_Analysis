clc        % Clear the command window
clear      % Clear all variables from the workspace
close all  % Close all figure windows

syms u(x,y) an(y)

% Define the function to be integrated
u(x,y) = sin(x) * sin(y);

% Compute the exact integral using symbolic integration
an(y) = int(u, x, [0 pi]);
ans = double(int(an, y, [0 pi]));

% Initialize function value matrix with ones
val = zeros(101) + 1;

syms x(i)
x(i) = i * (pi / 100);
y(i) = i * (pi / 100);

% Initialize Simpson's coefficient matrices
intmx = zeros(101);
for j = 0:100
    yc = j + 1;
    if j == 0 || j == 100
        intmx(:, yc) = intmx(:, yc) + 1;   % First and last terms
    elseif mod(j,2) == 1
        intmx(:, yc) = intmx(:, yc) + 4;   % Odd index terms
    else
        intmx(:, yc) = intmx(:, yc) + 2;   % Even index terms
    end
end

intr = zeros(101,1);
for i = 0:100
    xc = i + 1;
    if i == 0 || i == 100
        intr(xc,1) = intr(xc,1) + 1;
    elseif mod(i,2) == 1
        intr(xc,1) = intr(xc,1) + 4;
    else
        intr(xc,1) = intr(xc,1) + 2;
    end
end

% Compute function values at grid points
for i = 1:101
    val(:, i) = val(:, i) * sin(y(i-1));
    val(i, :) = val(i, :) * sin(x(i-1));
end

h = double(pi / 100);  % Step size

% Apply Simpson's rule for double integration
vh = val .* intmx;
vec = sum(vh, 2);
vec = (h / 3) * vec;
vec2 = intr .* vec;
ansh = sum(vec2);
ansh = (h / 3) * ansh;

% Uncomment to visualize function values
% mesh(val);
