
function [x_out, t_out, u_out] = FisherKolmogorov(nu, ro, K_cap, u0_func, L, Tf, Nx, Nt)
% FisherKolmogorov: Resuelve el modelo con condiciones de Neumann (ux=0)
% Entradas:
%   nu: Difusividad
%   ro: Tasa de crecimiento
%   K_cap: Capacidad de carga (Renombrado de 'k' para evitar error)
%   u0_func: Función anónima para la condición inicial
%   L, Tf: Dominio espacial [-L, L] y tiempo final
%   Nx, Nt: Número de nodos espacial y temporal

% --- 1. Datos y Adimensionalización ---
% Basado en las formulas del PDF
Adim = (ro * L^2) / nu;   % Parámetro adimensional lambda o A
T_adim = (Tf * nu) / (L^2); % Tiempo final adimensional

% --- 2. Discretización ---
x = linspace(-1, 1, Nx)';
t = linspace(0, T_adim, Nt)';

dx = 2 / (Nx - 1);        % Paso espacial adimensional (h)
dt = T_adim / (Nt - 1);   % Paso temporal adimensional

u = zeros(Nx, Nt);

% --- 3. Inicialización ---
% La condición inicial u0 depende de la variable dimensional (L*x)
% y se normaliza por K_cap para hacerla adimensional
u(:,1) = (1 / K_cap) * u0_func(L * x);

% --- 4. Matrices del Sistema (Esquema Semi-Implícito) ---

lamb = dt / (dx^2); % Número de Courant difusivo (lambda)

% Matriz de Difusión (Laplaciano implícito -1, 2, -1)
e = ones(Nx, 1);
% Diagonal inferior (-lamb), Principal (1 + 2*lamb), Superior (-lamb)
M_diff = spdiags([-lamb*e, (1 + 2*lamb)*e, -lamb*e], -1:1, Nx, Nx);

% --- 5. Iteraciones en tiempo ---
% IMPORTANTE: Usamos 'n' como índice, NO 'k'
for n = 1:(Nt - 1)
    u_curr = u(:, n); % Solución en el paso actual (n)

    % -- Parte de Reacción (Linealización Semi-Implícita) --
    % Término: Adim * u_curr * (1 - u_next)
    % Al pasar al LHS, sumamos a la diagonal: dt * Adim * u_curr

    vec_react = dt * Adim * u_curr;
    M_react = spdiags(vec_react, 0, Nx, Nx);

    % Matriz del sistema (LHS)
    A = M_diff + M_react;

    % Vector lado derecho (RHS)
    % Proviene de u^n/dt * dt = u^n más el término fuente linealizado
    b = u_curr .* (1 + dt * Adim);

    % --- 6. Condiciones de Contorno: Neumann Homogéneo (ux = 0) ---
    % Aproximación orden 1: u(1) = u(2) y u(N) = u(N-1)
    % Modificamos la matriz A y el vector b para forzar estas ecuaciones.

    % Nodo 1: u_1 - u_2 = 0
    A(1, :) = 0;        % Limpiar fila 1
    A(1, 1) = -1;
    A(1, 2) = 1;
    b(1) = 0;

    % Nodo N: u_N - u_{N-1} = 0
    A(Nx, :) = 0;       % Limpiar fila N
    A(Nx, Nx) = 1;
    A(Nx, Nx-1) = -1;
    b(Nx) = 0;

    % Resolución del sistema lineal
    u(:, n+1) = A \ b;
end

% --- 7. Recuperar variables dimensionales para la salida ---
x_out = L * x;
t_out = (L^2 / nu) * t;
u_out = K_cap * u;

end



