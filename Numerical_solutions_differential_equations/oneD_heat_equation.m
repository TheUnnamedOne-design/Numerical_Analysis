clc
clear all
syms xi(i) tj(j)

% Define spatial and time grid parameters
M=7;      % Number of spatial intervals
N=100;    % Number of time intervals
a=0;
b=1;
ti=0;
tf=1;

% Compute grid spacing
dx=(b-a)/M;
dt=(tf-ti)/N;

% Initialize solution matrix
mat1=zeros(N+1,M+1);

% Define spatial and time step locations
xi(i)=i*dx;
tj(j)=j*dt;

% Apply Dirichlet boundary conditions: u(0,t) = 0, u(1,t) = 0
for r=1:N+1
    mat1(r,M+1)=0;
    mat1(r,1)=0;
end

% Set initial condition: u(x,0) = sin(pi*x)
r=1;
for c=2:M+1
    mat1(r,c)=sin(pi*xi(c-1));
end

% Compute the stability parameter for explicit method
dr=dt/(dx^2);

% Apply explicit finite difference scheme
for r=2:N+1
    for c=2:M
        v1=mat1(r-1,c-1);
        v2=mat1(r-1,c);
        v3=mat1(r-1,c+1);
        mat1(r,c)=(dr*(v1+v3))+((1-2*dr)*(v2));
    end
end

% Compute analytical solution for comparison
mat2=zeros(N+1,M+1);
for r=1:N+1
    for c=1:M+1
        mat2(r,c)=(exp(-((pi^2)*tj(r-1))))*sin(pi*xi(c-1));
    end
end

% Extract second time step for visualization
mat3=zeros(1,M+1);
for i=1:M+1
    mat3(1,i)=mat1(2,i);
end

% Display numerical solution at the second time step
disp(mat3);
plot(mat3);

% Plot analytical solution
figure
xlabel("x")
ylabel("y")
contourf(transpose(mat2),200,'linecolor','non')
colormap(jet(256))
colorbar
title("Analytical Solution");
caxis([0,1])

% Plot numerical solution
figure
xlabel("x")
ylabel("y")
contourf(transpose(mat1),200,'linecolor','non')
colormap(jet(256))
colorbar
title("Numerical Solution");
caxis([0,1])
