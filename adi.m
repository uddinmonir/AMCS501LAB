% Parameters
N = 50; % Number of grid points in each direction
h = 1/(N+1); % Grid spacing
alpha = 1.0; % Parameter for ADI (can be tuned)

% Initialize grid
x = linspace(0, 1, N+2);
y = linspace(0, 1, N+2);
[X, Y] = meshgrid(x, y);

% Right-hand side function f(x, y)
f = @(x, y) sin(pi*x) .* sin(pi*y); % Example function

% Discretize f on the grid
F = f(X, Y);
F = F(2:end-1, 2:end-1); % Exclude boundary points
b = -h^2 * F(:); % Right-hand side vector

% Construct tridiagonal matrices for x and y directions
I = eye(N);
D = diag(-2*ones(N,1)) + diag(ones(N-1,1), 1) + diag(ones(N-1,1), -1);
A_x = kron(I, D); % Kronecker product for x-direction
A_y = kron(D, I); % Kronecker product for y-direction
I = eye(size(A_x));
% ADI iteration
u = zeros(N*N, 1); % Initial guess
max_iter = 1000; % Maximum number of iterations
tol = 1e-6; % Tolerance for convergence

for iter = 1:max_iter
    % Step 1: Solve in x-direction
    u = (I - alpha * A_x) \ (u + alpha * b);
    
    % Step 2: Solve in y-direction
    u = (I - alpha * A_y) \ (u + alpha * b);
    
    % Check for convergence
    residual = norm(A_x * u + A_y * u - b, 2);
    if residual < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end

% Reshape solution to 2D grid
U = reshape(u, N, N);

% Plot the solution
figure;
surf(X(2:end-1, 2:end-1), Y(2:end-1, 2:end-1), U);
title('Solution using ADI Method');
xlabel('x');
ylabel('y');
zlabel('u(x, y)');