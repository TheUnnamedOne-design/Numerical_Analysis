clc;
clear all;
close all;
tic;

% Define symbolic variables
syms x(i) y(j) t(n) u(Y)

% Define grid and time steps
p = 50;  % Number of x divisions
q = 50;  % Number of y divisions
td = 10^3;  % Number of time steps

% Define spatial and temporal boundaries
xi = 0; xf = 1;
yi = 0; yf = 1;
ti = 0; tf = 1;

% Compute grid step sizes
dx = (xf - xi) / p;
dy = (yf - yi) / q;
dt = (tf - ti) / td;

% Define coordinate arrays
x(i) = xi + i * dx;
y(j) = yi + j * dy;
t(n) = ti + n * dt;

% Define coefficients for finite difference scheme
ax = dt / (4 * dx);
by = dt / (2 * dy^2);
bx = dt / (2 * dx^2);

% Compute initial condition for u(Y)
t0 = 0;
k = 10;
yp = asinh(t0 * cosh(k)) / k;  % Compute yp
Pe = 10^3;  % Peclet number

% Define piecewise function for u(Y)
u(Y) = piecewise(Y < yp, (1 - (cosh(k * yp) / cosh(k)) - (k * (1 - yp)) * (sinh(k * yp) / cosh(k))), ...
                  Y >= yp, (1 - (cosh(k * Y) / cosh(k)) - (k * (1 - Y)) * (sinh(k * yp) / cosh(k))));

% Initialize solution grid
fs = zeros(p + 3, q + 3, td + 1);

% Apply Dirac delta function at x(i)
for i = 0:p
    xc = i + 2;
    val = dirac(x(i));
    if val == Inf
        val = 0;
    end
    fs(xc:xc, 2:q+2, 1) = val;
end

% Set initial condition
fs(:,:,1) = 1;

% Construct coefficient matrix for the linear system
pr = (p + 1) * (q + 1);
ls = zeros(pr);

for i = 1:pr
    ctrl = i - 1;
    yi = mod(ctrl, q + 1);
    yc = abs(y(yi));
    
    % Define coefficients based on grid location
    if mod(i, q + 1) == 0  % Boundary at y_max
        a = -bx - (Pe * u(yc) * ax);
        b = 2 * (-by);
        c = 1 + 2 * bx + 2 * by;
        d = 0;
        e = -bx + (Pe * u(yc) * ax);
    elseif mod(i, q + 1) == 1  % Boundary at y_min
        a = -bx - (Pe * u(yc) * ax);
        b = 0;
        c = 1 + 2 * bx + 2 * by;
        d = 2 * (-by);
        e = -bx + (Pe * u(yc) * ax);
    else  % Interior points
        a = -bx - (Pe * u(yc) * ax);
        b = -by;
        c = 1 + 2 * bx + 2 * by;
        d = -by;
        e = -bx + (Pe * u(yc) * ax);
    end

    % Modify boundary conditions
    if i <= q + 1
        e = a + e;
        a = 0;
    end
    if i >= pr - q - 1
        a = a + e;
        e = 0;
    end

    % Assign values to coefficient matrix
    ls(i, i) = c;
    if i - 1 >= 1
        ls(i, i - 1) = b;
    end
    if i + 1 <= pr
        ls(i, i + 1) = d;
    end
    if i - q - 1 >= 1
        ls(i, i - q - 1) = a;
    end
    if i + q + 1 <= pr
        ls(i, i + q + 1) = e;
    end
end

% Initialize vectors for finite difference computations
rs = zeros(pr, 1);
al = ones(p+1, q+1) * bx;
ga = ones(p+1, q+1) * by;
ga(:, q+1) = ga(:, q+1) + by;
oh = ones(p+1, q+1) * (1 - 2 * bx - 2 * by);
si = ones(p+1, q+1) * by;
si(:, 1) = si(:, 1) + by;
la = ones(p+1, q+1) * bx;

% Compute correction terms based on boundary conditions
for j = 0:q
    yc = j + 1;
    Yh = abs(y(j));
    al(:, yc) = al(:, yc) + Pe * u(Yh) * ax;
    la(:, yc) = la(:, yc) - Pe * u(Yh) * ax;
end

% Time-stepping loop
for n = 1:td
    tc = n + 1;

    % Compute right-hand side for time step
    d1 = double(al .* fs(1:p+1, 2:q+2, tc-1));
    d2 = double(ga .* fs(2:p+2, 1:q+1, tc-1));
    d3 = double(oh .* fs(2:p+2, 2:q+2, tc-1));
    d4 = double(si .* fs(2:p+2, 3:q+3, tc-1));
    d5 = double(la .* fs(3:p+3, 2:q+2, tc-1));
    d0 = d1 + d2 + d3 + d4 + d5;

    % Convert to column vector for solving
    d0 = double(transpose(d0));
    rs = d0(:);

    % Solve system of equations
    soln = double(inv(ls) * rs);
    holder = reshape(soln, [p+1, q+1]);
    holder = transpose(holder);
    fs(2:p+2, 2:q+2, tc) = holder;

    % Store a representative value for plotting
    Cholder(1, n) = fs(7, 7, n);
end

toc;

% Compute derived quantities for visualization
dzt = xzC ./ xC;
dz2t = xz2C ./ xC;
dzt(1,1) = 0;
dz2t(1,1) = 0;

tv = ti + dt : dt : tf;
dztsq = dzt .* dzt;
vhold = dz2t - dztsq;
dezt = vhold(1:1, 2:end) - vhold(1:1, 1:td);

% Plot results
ax2 = nexttile;
semilogx(tv, Cholder, 'b-', 'LineWidth', 2);
hold on;
