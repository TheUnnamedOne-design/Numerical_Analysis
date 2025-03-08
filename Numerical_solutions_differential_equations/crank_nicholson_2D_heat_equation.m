

% This code contains the method of solution and 
% the plotting of said solution in matlab. The two dimensional heat equation has been solved in the following code
%Refer to Problem in pdf 2D_heat_equation to refer to the problem definition


clc  % Clear command window
clear all  % Clear all variables
close all  % Close all figures

syms x(i) y(j) t(n)  % Define symbolic variables for spatial and time grids

% Grid parameters
p = 10;  % Number of spatial divisions in x-direction
q = 10;  % Number of spatial divisions in y-direction
r = 30;  % Number of time steps
alphac = 0.05;  % Thermal diffusivity

% Domain boundaries
ti = 0; tf = 1;  % Time domain
xi = 0; xf = 1;  % X domain
yi = 0; yf = 1;  % Y domain

% Step sizes
dx = (xf - xi) / p;  
dy = (yf - yi) / q;
dt = (tf - ti) / r;

% Discretized space and time
x(i) = xi + i * dx;
y(j) = yi + j * dy;
t(n) = ti + n * dt;

% Crank-Nicholson coefficients
rx = (alphac * dt) / (dx)^2;
ry = (alphac * dt) / (dy)^2;
al = -rx / 2;
bl = -ry / 2;
cl = (1 + rx + ry);

% Initialize coefficient matrix (left-hand side)
ls = zeros(((p+1)*(q+1)), ((p+1)*(q+1)));  
rs = zeros(((p+1)*(q+1)), 1);  % Right-hand side matrix

% Constructing the coefficient matrix
for i = 1:(p+1)*(q+1)
    if(i+1 <= (p+1)*(q+1))
        ls(i, i+1) = bl;  % Right neighbor coefficient
    end
    if(i-1 >= 1)
        ls(i, i-1) = bl;  % Left neighbor coefficient
    end
    if(i+q+1 <= (p+1)*(q+1))
        ls(i, i+q+1) = al;  % Top neighbor coefficient
    end
    if(i-q-1 >= 1)
        ls(i, i-q-1) = al;  % Bottom neighbor coefficient
    end
    ls(i, i) = cl;  % Central coefficient

    % Apply boundary conditions to ensure correct indexing
    if(mod(i, q+1) == 0 && i ~= (p+1)*(q+1))
        ls(i, i+1) = 0;
    end
    if(mod(i, q+1) == 1 && i ~= 1)
        ls(i, i-1) = 0;
    end
end

% Initialize solution matrix for temperature field
fs = zeros(p+1, q+1, r+1);

% Apply initial boundary conditions
for i = 1:p+1
    for j = 1:q+1
        if(i == 1 || i == p+1)
             fs(i, j, 1) = 1;  % Dirichlet boundary condition (fixed temperature)
        end

    end
end

% Time-stepping loop
for n = 1:r
    tc = n+1;  % Current time step index
    count = 1; % Counter for coefficient matrix index
    
    % Construct right-hand side matrix (explicit part)
    for i = 0:p
        xc = i+1;
        for j = 0:q
            yc = j+1;
            d0 = (1 - rx - ry) * fs(xc, yc, tc-1);
            d1 = 0; d2 = 0; d3 = 0; d4 = 0;
            
            % Compute contributions from neighboring points
            if(j > 0) 
                d1 = (ry / 2) * fs(xc, yc-1, tc-1);
            end
            if(j < q)
                d2 = (ry / 2) * fs(xc, yc+1, tc-1);
            end
            if(i > 0)
                d3 = (rx / 2) * fs(xc-1, yc, tc-1);
            end
            if(i < p)
                d4 = (rx / 2) * fs(xc+1, yc, tc-1);
            end

            % Right-hand side matrix for Crank-Nicholson method
            dn = d0 + d1 + d2 + d3 + d4;
            rs(count, 1) = dn;
            count = count + 1;
        end
    end

    % Solve the linear system
    soln = inv(ls) * rs;

    % Store solution into fs matrix
    count = 1;
    for i = 0:p
        xc = i+1;
        for j = 0:q
            yc = j+1;
            fs(xc, yc, tc) = soln(count, 1);
            if(i == 0 || i == p)  % Enforce boundary condition at every time step
                fs(xc, yc, tc) = 1;
            end
            count = count + 1;
        end
    end
end

fs;  % Final solution matrix

% Visualization of temperature evolution
figure
xlabel("x")
ylabel("y")
for i = 1:r+1
    contourf(fs(:, :, i), 200, 'linecolor', 'non')  % Contour plot
    colormap(jet(256))
    colorbar
    val = double(t(i-1));
    vh = round(val, 3);
    
    % Generate dynamic title based on time step
    str = "The plot at time ";
    str2 = compose('%.3f', val);
    str3 = "s";
    sh = strcat(str, str2, str3);
    title(sh);
    
    caxis([0,1])  % Set color scale
    
    % Pause to visualize time evolution
    if(i == 1)
        pause(5);
    else
        pause(dt);
    end
    hold on;
end
