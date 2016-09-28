close all; clear; clc;
w = warning ('off','all');
%format long;
% PROBLEMA DE FLUJO ÓPTIMO

% min cX
%   s.t. Nx = b
%        tx <= T
%        x >= 0

% initializo la matriz nodo/arco
% N = [1  1  0  0  0  0  0  0  0  0; ...
%     -1  0  1  1  0  0  0  0  0  0; ...
%      0 -1  0  0  1  1  1  0  0  0; ...
%      0  0 -1  0  0 -1  0  1  1  0; ...
%      0  0  0  0  -1  0 -1 -1  0  1; ...
%      0  0  0 -1  0  0  0  0 -1 -1];
 
 N = [ 1  1  0  0  0  0  0;
      -1  0  1  1  0  0  0;
       0 -1  0  0  1  0  0;
       0  0 -1  0  0  1  0;
       0  0  0  0 -1  0  1;
       0  0  0 -1  0 -1 -1;]
% b: flujo que entra y flujo que sale
b = [1 0 0 0 0 -1];
% c: costos
%c = [1 10 1 2 1 5 12 10 1 2];
c = [2 1 2 5 2 1 2];
% t: tiempos
t = [3 1 3 1 3 3 5];
% t: cota de tiempo
T = 8;

% Con esta sentencia me aseguro de usar el método simplex, que como recorre
% los vértices, en caso de que haya 2 caminos óptimos iguales, va a
% devolvernos uno solo. Si no usara esto, obtendría flujos fraccionarios.
options = optimoptions('linprog', 'Algorithm', 'Simplex', 'Display', 'off');

% Resuelvo por cutting planes

lambdas = [];           %
f_lambdas = [];         %
f_lambdas_real = [];    % estructuras para acumular valores y graficar
xis = [];               %
x_stars = [];           %

lambda_0 = 1;
lambda_i = lambda_0;

epsilon = 1e-4;

% precalculo los valores de x_star, xi y f_lambda para inicializar
lambdas = [lambda_0 lambdas];

[x_star_i, f_lambda_real_i] = linprog(c + lambda_0 * t, [], [], N, b, zeros(1,7), [],[],options); 
f_lambda_i = (c + lambda_0 * t) * x_star_i - lambda_0 * T;
f_lambdas = [f_lambdas f_lambda_i];
f_lambdas_real = [f_lambdas_real f_lambda_real_i];
xi_i = t * x_star_i - T;
xis = [xis xi_i];
x_stars = [x_stars x_star_i];

A = [-1 -xi_i]; 
B = [f_lambda_i - xi_i * lambda_i];

iter = 0;

while (abs(f_lambda_real_i - f_lambda_i) > epsilon && iter < 15)
    
    % busco mi nuevo lambda minimizando %CUTTING PLANES
    o = linprog([1 0], A, B, [], [], [-20 -20], [20 20], [], options);
    lambda_i = o(end);  
    
    lambda_i_prev = lambda_i;
    % MODELO ORIGINAL
    [x_star_i, f_lambda_real_i] = linprog(c + lambda_i * t, [], [], N, b, zeros(1,7),[],[],options);
    xi_i = t * x_star_i - T;
    t * x_star_i
    f_lambda_i = (c + lambda_i * t) * x_star_i - lambda_i * T;
    
    f_lambdas = [f_lambdas f_lambda_i];                   %
    f_lambdas_real = [f_lambdas_real f_lambda_real_i ];   %
    xis = [xis xi_i];                                     % estruct. para 
    x_stars = [x_stars; x_star_i];                        % acumular
    lambdas = [lambdas lambda_i];                         %
    
    A = cat(1, A, [-1 -xi_i]);
    B = cat(2, B, f_lambda_i - xi_i * lambda_i);
    
    iter = iter + 1;
end

figure, plot(f_lambdas); title('f(lambda) estimado');
figure, plot(f_lambdas_real); title('f(lambda) real');
figure, plot(xis); title('xis');
figure, plot(x_star_i); title('x*');
figure, plot(lambdas); title('lambdas');
