close all, clear, clc;
w = warning ('off','all');

% PROBLEMA DE FLUJO ÓPTIMO

% min cX
%   s.t. Nx = b
%        tx <= T
%        x >= 0

% initializo la matriz nodo/arco
N = [1 1 0 0 0 0 0 0 0 0; ...
    -1 0 1 1 0 0 0 0 0 0; ...
    0 -1 0 0 1 1 1 0 0 0; ...
    0 0 -1 0 0 -1 0 1 1 0; ...
    0 0 0 -1 0 0 -1 -1 0 1; ...
    0 0 0 0 0 0 0 0 -1 -1];

% b: flujo que entra y flujo que sale
b = [1 0 0 0 0 -1];
% c: costos
c = [1 10 1 2 1 5 12 10 1 2];
% t: tiempos
t = [10 3 1 3 2 7 3 1 7 2];
% t: cota de tiempo
T = 10;

% Con esta sentencia me aseguro de usar el método simplex, que como recorre
% los vértices, en caso de que haya 2 caminos óptimos iguales, va a
% devolvernos uno solo. Si no usara esto, obtendría flujos fraccionarios.
options = optimoptions('linprog', 'Algorithm', 'Simplex', 'Display', 'off');

% Resuelvo por cutting planes

lambda_0 = 0;
lambda_i_prev = -Inf;
lambda_i = lambda_0;

epsilon = 1e-4;

% precalculo los valores de x_star, xi y f_lambda para inicializar
lambdas = lambda_i;
[x_stars, xis, f_lambdas] = oraculo(lambda_i, c, t, T, N, b);

% itero hasta que la diferencia entre los lambdas sea menor a epsilon
while abs(lambda_i - lambda_i_prev) > epsilon
    
    % guardo el valor actual de lambda como valor anterior
    lambda_i_prev = lambda_i;
    
    % llamo al oráculo para que me de un nuevo valor de x_star y de
    % subgradiente xi
    [x_star_i, xi_i, f_lambda_i] = oraculo(lambda_i, c, t, T, N, b);

    % concateno los valores del nuevo plano
    x_stars = cat(2, x_stars, x_star_i);
    xis = cat(1, xis, xi_i);
    f_lambdas = cat(1, f_lambdas, f_lambda_i);
    
    % busco mi nuevo lambda minimizando
    A = cat(2, -ones(length(xis), 1), -xis);
    b = cat(1, f_lambdas - xis .* lambdas);
    output = linprog([1 0], A, b, [], [], -Inf, 0);
    
    lambda_i = output(end);
    
end






% %CÓDIGO VIEJO
% % Resuelvo el problema sin restricciones de tiempo
% disp('Problema sin restricciones de tiempo:');
% xsol = linprog(c, [], [], N, b, zeros(size(N,1),1), [], [], options);
% 
% % defino phi de lambda
% 
% lambda_vals = 0:0.1:10;
% phi_lambda = zeros(length(lambda_vals));
% for i=1:length(lambda_vals)
%     lambda=lambda_vals(i);
%     [sol, fval] = linprog((c + lambda * t'), [], [], N, b, zeros(size(N,1),1), [], [], options);
%     phi_lambda(i) = -lambda * T + fval;
% end
% 
% figure, plot(lambda_vals, phi_lambda);
% title('Función objetivo del problema dual');
% xlabel('\lambda values');
% ylabel('\phi(\lambda)');
% box on
% grid on

    