close all;
clc;

% Physical constants
mu_0 = 4 * pi * 1e-7;
eps_0 = 8.854e-12;
eps_r = 2.3;
mu_r = 1;

% Grid parameters
Nx = 101;
Ny = Nx;
Lx = 3; % mm
Ly = Lx;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dx_metric = dx * 1e-3;
dy_metric = dy * 1e-3;

V = zeros(Nx, Ny);
Jz = zeros(Nx, Ny);

x = linspace(-Lx / 2, Lx / 2, Nx);
y = linspace(-Ly / 2, Ly / 2, Ny);
[X, Y] = meshgrid(x, y);

% Elliptical conductor dimensions in mm
a_in = 0.5;
b_in = 0.25;
a_out = 1.2;
b_out = 1;
a_mid = (a_out + a_in) / 2;
b_mid = (b_out + b_in) / 2;
x_offset = 0.00;

inside_inner = ((X - x_offset) / a_in).^2 + (Y / b_in).^2 <= 1;
outside_outer = (X / (a_out + a_in / 2)).^2 + ...
    (Y / (b_out + b_in / 2)).^2 <= 1;
inside_outer_boundary = (X / a_out).^2 + (Y / b_out).^2 >= 1;
inside_outer = outside_outer & inside_outer_boundary;

inside_dielectric_inner = ((X - x_offset) / a_mid).^2 + ...
    (Y / b_mid).^2 <= 1;
inside_dielectric_outer = ((X - x_offset) / a_mid * 1.1).^2 + ...
    (Y / b_mid * 1.1).^2 >= 1;
inside_dielectric_layer = inside_dielectric_inner & inside_dielectric_outer;

% Conductor excitation
Voltage = 2;
Current = 2;
V_in = Voltage;
V_out = 0;
J_in = Current / (sum(inside_inner(:)) * dx_metric * dy_metric);
J_out = -Current / (sum(inside_outer(:)) * dx_metric * dy_metric);

V(inside_inner) = V_in;
V(inside_outer) = V_out;
Jz(inside_inner) = J_in;
Jz(inside_outer) = J_out;

% Solve for electric potential and magnetic vector potential
Phi = FDM_Solver(Nx, V, inside_inner, inside_outer, 'E');
Az = FDM_Solver(Nx, mu_0 * mu_r * Jz * dx_metric * dx_metric, ...
    inside_inner, inside_outer, 'B');

% Electric field
[Ex, Ey] = gradient(-Phi, dx_metric, dy_metric);

% Magnetic flux density from the z-directed vector potential
[dAz_dx, dAz_dy] = gradient(Az, dx_metric, dy_metric);
Bx = dAz_dy;
By = -dAz_dx;

% 2D mesh and conductor locations
figure;
plot(X, Y, 'k.');
hold on;
plot(X', Y', 'k.');
plot(X(inside_inner), Y(inside_inner), 'ro', ...
    'MarkerSize', 1, 'MarkerFaceColor', 'r');
plot(X(inside_outer), Y(inside_outer), 'ro', ...
    'MarkerSize', 1, 'MarkerFaceColor', 'b');
title('2D Mesh mit Leitern');
xlabel('x [mm]');
ylabel('y [mm]');
axis equal;
grid on;
x_margin = (max(x) - min(x)) * 0.05;
y_margin = (max(y) - min(y)) * 0.05;
xlim([min(x) - x_margin, max(x) + x_margin]);
ylim([min(y) - y_margin, max(y) + y_margin]);

% Mesh with the integration contour for electric flux
figure;
plot(X, Y, 'k.');
hold on;
plot(X', Y', 'k.');
plot(X(inside_inner), Y(inside_inner), 'ro', ...
    'MarkerSize', 1, 'MarkerFaceColor', 'r');
plot(X(inside_outer), Y(inside_outer), 'ro', ...
    'MarkerSize', 1, 'MarkerFaceColor', 'b');
plot(X(inside_dielectric_layer), Y(inside_dielectric_layer), 'bo', ...
    'MarkerSize', 1, 'MarkerFaceColor', 'b');
title('Huellenintegral');
xlabel('x [mm]');
ylabel('y [mm]');
axis equal;
grid on;
x_margin = (max(x) - min(x)) * 0.05;
y_margin = (max(y) - min(y)) * 0.05;
xlim([min(x) - x_margin, max(x) + x_margin]);
ylim([min(y) - y_margin, max(y) + y_margin]);

% Potential and field plots
figure;
contourf(X, Y, Phi, 50);
colorbar;
title('Potentialverteilung [V]');
xlabel('x');
ylabel('y');
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

figure;
quiver(X, Y, Ex, Ey);
title('Elektrisches Feld [V/m]');
xlabel('x');
ylabel('y');
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

figure;
quiver(X, Y, Bx, By);
title('Magnetische Flussdichte [T]');
xlabel('x');
ylabel('y');
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

figure;
contourf(X, Y, Jz, 50);
colorbar;
title('Eingepragte Stromdichte [A/m^2]');
xlabel('x');
ylabel('y');
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

V(inside_inner) = 1;
V(inside_outer) = -1;
figure;
contourf(X, Y, V, 50);
colorbar;
title('Eingepragtes Potential [V]');
xlabel('x');
ylabel('y');
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

% Extract one-dimensional field profiles at y = 0
m = size(By, 1);
mid_row = round(m / 2);
B_1D = By(mid_row, :);
E_1D = Ex(mid_row, :);

figure;
ax1 = axes();
plot(ax1, x, B_1D, 'r-x', 'LineWidth', 1.5);
xlabel('x');
ylabel('Magnetische Flussdichte By (T)');
title('Magnetische Flussdichte und elektrisches Feld bei Y = 0');
grid on;
hold on;
plot(ax1, [a_in + x_offset, a_in + x_offset], ylim(ax1), ...
    'k-', 'LineWidth', 1.5);
plot(ax1, [-(a_in - x_offset), -(a_in - x_offset)], ylim(ax1), ...
    'k-', 'LineWidth', 1.5);
plot(ax1, [a_out, a_out], ylim(ax1), 'k-', 'LineWidth', 1.5);
plot(ax1, [-a_out, -a_out], ylim(ax1), 'k-', 'LineWidth', 1.5);

ax2 = axes();
plot(ax2, x, E_1D, 'b-o', 'LineWidth', 1.5);
ylabel('Elektrisches Feld Ex (V/m)');
set(ax2, 'Color', 'none', 'YAxisLocation', 'right', 'XColor', 'none');
legend(ax1, 'B (T)', 'Location', 'NorthWest');
legend(ax2, 'E (V/m)', 'Location', 'NorthEast');
hold off;

% Capacitance per unit length from electric flux
sigma = eps_0 * eps_r * sqrt(Ex.^2 * dx_metric^2 + Ey.^2 * dy_metric^2);
Q = 0;
for i = 1:Nx
    for j = 1:Ny
        if inside_dielectric_layer(i, j)
            if i > Nx / 2
                if ~inside_dielectric_layer(i + 1, j)
                    Q = Q + sigma(i, j);
                end
            elseif ~inside_dielectric_layer(i - 1, j)
                Q = Q + sigma(i, j);
            end

            if j > Ny / 2
                if ~inside_dielectric_layer(i, j + 1)
                    Q = Q + sigma(i, j);
                end
            elseif ~inside_dielectric_layer(i, j - 1)
                Q = Q + sigma(i, j);
            end
        end
    end
end
U = V_in - V_out;
C_prime = Q / U;

% Inductance per unit length from magnetic flux
[~, idx_offset] = min(abs(x - x_offset));
B_all = B_1D(idx_offset:end);
Mag_Phi = sum(B_all(:)) * dx_metric;
L_prime = Mag_Phi / Current;
Z = sqrt(L_prime / C_prime);

disp(['Capacitance per unit length C'' = ', num2str(C_prime), ' F/m']);
disp(['Inductance per unit length L'' = ', num2str(L_prime), ' H/m']);
disp(['Impedance Z = ', num2str(Z), ' Ohm']);

% Circular-coax reference values for comparison
r_in = a_in * 1e-3;
r_out = a_out * 1e-3;
C_prime_theory = (2 * pi * eps_0 * eps_r) / log(r_out / r_in);
L_prime_theory = (mu_0 * mu_r) / (2 * pi) * log(r_out / r_in);
Z_theory = sqrt(L_prime_theory / C_prime_theory);

disp(['Reference capacitance C'' = ', ...
    num2str(C_prime_theory * 1e9), ' nF/m']);
disp(['Reference inductance L'' = ', ...
    num2str(L_prime_theory * 1e6), ' uH/m']);
disp(['Reference impedance Z = ', num2str(Z_theory), ' Ohm']);

% Summary figure
figure;
subplot(3, 2, 1);
contourf(X, Y, Phi, 20);
colorbar;
title('Potentialverteilung [V]');
xlabel('x');
ylabel('y');
axis equal;

subplot(3, 2, 2);
contourf(X, Y, Jz, 50);
colorbar;
title('Stromdichte [A/m^2]');
xlabel('x');
ylabel('y');
axis equal;

step = 3;
X_sub = X(1:step:end, 1:step:end);
Y_sub = Y(1:step:end, 1:step:end);
Ex_sub = Ex(1:step:end, 1:step:end);
Ey_sub = Ey(1:step:end, 1:step:end);

subplot(3, 2, 3);
quiver(X_sub, Y_sub, Ex_sub, Ey_sub, 'AutoScale', 'on', ...
    'AutoScaleFactor', 2, 'LineWidth', 1.5);
title('Elektrisches Feld [V/m]');
xlabel('x');
ylabel('y');
axis equal;

Bx_sub = Bx(1:step:end, 1:step:end);
By_sub = By(1:step:end, 1:step:end);
subplot(3, 2, 4);
quiver(X_sub, Y_sub, Bx_sub, By_sub, 'AutoScale', 'on', ...
    'AutoScaleFactor', 2, 'LineWidth', 1.5);
title('Magnetische Flussdichte [T]');
xlabel('x');
ylabel('y');
axis equal;

subplot(3, 2, 5);
plot(x, B_1D, 'r-x', 'LineWidth', 1.5);
xlabel('x');
ylabel('Magnetische Flussdichte By (T)');
title('Magnetische Flussdichte bei Y = 0');
grid on;
hold on;
plot([a_in + x_offset, a_in + x_offset], ylim, 'k-', 'LineWidth', 1.5);
plot([-(a_in - x_offset), -(a_in - x_offset)], ylim, ...
    'k-', 'LineWidth', 1.5);
plot([a_out, a_out], ylim, 'k-', 'LineWidth', 1.5);
plot([-a_out, -a_out], ylim, 'k-', 'LineWidth', 1.5);

subplot(3, 2, 6);
plot(x, E_1D, 'b-o', 'LineWidth', 1.5);
ylabel('Elektrisches Feld Ex (V/m)');
xlabel('x');
title('Elektrisches Feld bei Y = 0');
grid on;
hold on;
plot([a_in + x_offset, a_in + x_offset], ylim, 'k-', 'LineWidth', 1.5);
plot([-(a_in - x_offset), -(a_in - x_offset)], ylim, ...
    'k-', 'LineWidth', 1.5);
plot([a_out, a_out], ylim, 'k-', 'LineWidth', 1.5);
plot([-a_out, -a_out], ylim, 'k-', 'LineWidth', 1.5);

function Phi_mat = FDM_Solver(n, source, inside_inner, inside_outer, field)
    N = n * n;
    A = sparse(N, N);
    b = zeros(N, 1);

    for i = 2:n - 1
        for j = 2:n - 1
            idx = (j - 1) * n + i;
            if (~inside_inner(i, j) && ~inside_outer(i, j)) || field == 'B'
                A(idx, idx) = 4;
                A(idx, idx - 1) = -1;
                A(idx, idx + 1) = -1;
                A(idx, idx - n) = -1;
                A(idx, idx + n) = -1;
            else
                A(idx, idx) = 1;
            end
            b(idx) = source(i, j);
        end
    end

    % Neumann boundary conditions at the domain edges
    for i = 1:n
        for j = 1:n
            idx = (j - 1) * n + i;
            if i == 1
                A(idx, idx) = 1;
                A(idx, idx + 1) = 1;
                b(idx) = 0;
            elseif i == n
                A(idx, idx) = 1;
                A(idx, idx - 1) = 1;
                b(idx) = 0;
            elseif j == 1
                A(idx, idx) = 1;
                A(idx, idx + n) = 1;
                b(idx) = 0;
            elseif j == n
                A(idx, idx) = 1;
                A(idx, idx - n) = 1;
                b(idx) = 0;
            end
        end
    end

    tol = 1e-10;
    maxit = 1500;
    [Phi_vec, ~, ~] = gmres(A, b, [], tol, maxit);
    Phi_mat = reshape(Phi_vec, n, n);
end