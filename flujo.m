close all, clear, clc;
w = warning ('off','all');

% PROBLEMA DE FLUJO ÓPTIMO

% min cX
%   s.t. Nx = b
%        x >= 0

% initializo la matriz nodo/arco
N = [1 1 0 0 0 0 0; ...
    -1 0 1 1 0 0 0; ...
    0 -1 0 0 1 0 0; ...
    0 0 -1 0 0 1 0; ...
    0 0 0 0 -1 0 1; ...
    0 0 0 -1 0 -1 -1];

% b: flujo que entra y flujo que sale
b = [1 0 0 0 0 -1];
% c: costos
c = [1 1 2 6 2 2 3]';
% t: tiempos
t = [1 1 2 2 1 1 2];

T = 3.5;

% Con esta sentencia me aseguro de usar el método simplex, que como recorre
% los vértices, en caso de que haya 2 caminos óptimos iguales, va a
% devolvernos uno solo. Si no usara esto, obtendría flujos fraccionarios.
options = optimoptions('linprog', 'Algorithm', 'Simplex', 'Display', 'off');

% Resuelvo el problema sin restricciones de tiempo
disp('Problema sin restricciones de tiempo:');
xsol = linprog(c, [], [], N, b, zeros(size(N,1),1), [], [], options)

% Resuelvo el problema con restricciones de tiempo
disp('Problema con restricciones de tiempo:');
xsol_con_tiempo = linprog(c, t, T, N, b, zeros(size(N,1),1), [], [], options)

% Grafico los valores de la función objetivo \phi(\lambda) del problema
% dualizado
disp('Grafico los valores de mi función objetivo');
lambda_vals = 0:0.1:10;
phi_lambda = zeros(size(lambda_vals));
for i = 1 : length(lambda_vals)
    lambda = lambda_vals(i);
    [sol, fval] = linprog((c + lambda * t'), [], [], N, b, zeros(size(N,1),1), [], [], options);
    phi_lambda(i) = -lambda * T + fval;
end
figure, plot(lambda_vals, phi_lambda);
title('Función objetivo del problema dual');
xlabel('\lambda values');
ylabel('\phi(\lambda)');
box on
grid on

% Minimizamos utilizando gradient descent
disp('Obtenemos el máximo usando gradient descent');
epsilon = 0.001;
fval_(1) = 0;
lambda_raiz(1) = 0;
lambda_lineal(1) = 0;
lambda_cuadrado(1) = 0;
i = 2;
converge = false;
while ((~converge) || (i < 10)) && (i<300)
    lambda_raiz = cat(1, lambda_raiz, lambda_raiz(i-1) + (1/sqrt(i-1)) * (-T + t * linprog((c + lambda_raiz(i-1) * t'), [], [], N, b, zeros(size(N,1),1), [], [], options)));
    lambda_lineal = cat(1, lambda_lineal, lambda_lineal(i-1) + (1/(i-1)) * (-T + t * linprog((c + lambda_lineal(i-1) * t'), [], [], N, b, zeros(size(N,1),1), [], [], options)));
    lambda_cuadrado = cat(1, lambda_cuadrado, lambda_cuadrado(i-1) + (1/(i-1)^2) * (-T + t * linprog((c + lambda_cuadrado(i-1) * t'), [], [], N, b, zeros(size(N,1),1), [], [], options)));
    converge = abs(lambda_raiz(i) - lambda_raiz(i-1)) < epsilon; 
    i = i + 1;
end
hold on
plot(lambda_raiz, zeros(size(lambda_raiz)), '*');
legend({'Función objetivo', '\lambda_n'});

% Graficamos la evolución del error para 2 tamaños de paso \alpha distintos
disp('Graficamos la evolución del error para 2 tamaños de paso \alpha distintos');
lambda_raiz_error = ones(size(lambda_raiz)) * 2 - lambda_raiz;
lambda_lineal_error = ones(size(lambda_lineal)) * 2 - lambda_lineal;
lambda_cuadrado_error = ones(size(lambda_cuadrado)) * 2 - lambda_cuadrado;
figure, plot(1:length(lambda_raiz_error), lambda_raiz_error);
hold on
plot(1:length(lambda_lineal_error), lambda_lineal_error);
hold on
plot(1:length(lambda_cuadrado_error), lambda_cuadrado_error);
xlabel('Iteración');
ylabel('Diferencia respecto a 2');
title('Evolución del error en \lambda_n');
legend({'\alpha = 1/sqrt(n)', '\alpha = 1/n', '\alpha = 1/n^2'}, 'Location', 'northeast'); 
box on
grid on;
