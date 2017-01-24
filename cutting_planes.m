close all; clear; clc;
w = warning ('off','all');
%format long;
% PROBLEMA DE FLUJO ?PTIMO

% min cX
%   s.t. Nx = b
%        tx <= T
%        x >= 0

keys = {'1-2', '1-3', '2-4', '2-6', '3-5', '4-6', '5-6'};

% initializo la matriz nodo/arco
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

lower_bounds = [-Inf 0];






% Con esta sentencia me aseguro de usar el m?todo simplex, que como recorre
% los v?rtices, en caso de que haya 2 caminos ?ptimos iguales, va a
% devolvernos uno solo. Si no usara esto, obtendr?a flujos fraccionarios.
options = optimoptions('linprog', 'Algorithm', 'Simplex', 'Display', 'off');

% Resuelvo por cutting planes

lambdas = [];           %
f_lambdas = [];         %
phi_lambdas = [];    % estructuras para acumular valores y graficar
xis = [];               %
x_stars = [];           %

epsilon = 1e-4;

%%

lambda_0 = 0;%rand;
lambda_i = lambda_0;

% precalculo los valores de x_star, xi y f_lambda para inicializar
lambdas = [lambda_0 lambdas];

% resuelvo:
%
% min c' x + \lambda * (tx - T)
% x>=0
%
[x_star_i, phi_lambda_i_parcial] = linprog(c + lambda_0 * t, [], [], N, b, zeros(1,size(N,2)), [],[],options); 
phi_lambda_i = phi_lambda_i_parcial - lambda_0 * T;

% f_lambda_i es el resultado de resolver mi cutting planes. Como todav?a no
% lo corr?, inicialmente lo seteo en un valor alto
f_lambda_i = 30000;

% Concateno todo para ir llevando un hist?rico
f_lambdas = [f_lambdas f_lambda_i];
phi_lambdas = [phi_lambdas phi_lambda_i];
xi_i = t * x_star_i - T;
xis = [xis xi_i];
x_stars = [x_stars x_star_i];

% Acumulo los subgradientes
A = [-1 -xi_i]; 
% Acumulo lo que necesito para las desigualdades que resuelvo por cutting
% planes
B = [phi_lambda_i - xi_i * lambda_i];

%%

iter = 0;

while ~(abs(phi_lambda_i - f_lambda_i) < epsilon)% && (iter < 30)
    
    
    %CUTTING PLANES
    [p, f_lambda_i] = linprog([1 0], A, B, [], [], lower_bounds, [], [], options);
    f_lambda_i = f_lambda_i * -1;
    lambda_i = p(end);  
    lambda_i_prev = lambda_i;
    
    % MODELO ORIGINAL
    [x_star_i, phi_lambda_i_parcial] = linprog(c + lambda_i * t, [], [], N, b, zeros(1,size(N,2)),[],[],options);
    phi_lambda_i = phi_lambda_i_parcial - lambda_i * T;
    
    xi_i = t * x_star_i - T;
    %f_lambda_i = (c + lambda_i * t) * x_star_i - lambda_i * T;
    
    f_lambdas = [f_lambdas f_lambda_i];                   %
    phi_lambdas = [phi_lambdas phi_lambda_i ];   %
    xis = [xis xi_i];                                     % estruct. para 
    x_stars = [x_stars; x_star_i];                        % acumular
    lambdas = [lambdas lambda_i];                         %
    
    A = cat(1, A, [-1 -xi_i]);
    B = cat(2, B, phi_lambda_i - xi_i * lambda_i);
    
    [ phi_lambda_i,   f_lambda_i ]
    
    iter = iter + 1;
end

if abs(phi_lambda_i - f_lambda_i) <= epsilon
    keys(find(x_star_i))
else
   disp('No hay solucion');
end

% figure
% 
% subplot(2,2,1);
% plot(f_lambdas, 'LineWidth', 2); 
% hold on
% plot(f_lambdas_real, 'LineWidth', 2)
% legend('f(\lambda) estimado', 'f(\lambda) real', 'location', 'southeast');
% xlabel('iteraci?n')
% ylabel('f');
% xlim([1 iter]);
% grid on
% title('Evoluci?n del valor de f(\lambda) estimado y real por iteraci?n');
% 
% subplot(2,2,3);
% plot(xis, 'LineWidth', 2); 
% xlim([1 iter]);
% grid on
% title('Evoluci?n de \xi por iteraci?n');
% xlabel('iteraci?n')
% ylabel('\xi')
% 
% subplot(2,2,4);
% bar(x_star_i); 
% title('Soluci?n');
% xlabel('Columna en matriz nodo/arco')
% ylabel('x^*');
% 
% subplot(2,2,2);
% plot(lambdas, 'LineWidth', 2);
% xlim([1 iter]);
% grid on
% title('Evoluci?n de \lambda por iteraci?n');
% xlabel('iteraci?n')
% ylabel('\lambda');


min_lambda = 0;
prec = 0.1;
max_lambda = 10;

dom_lambda = min_lambda : prec : max_lambda;
phi_lambda = zeros(length(dom_lambda), 1);
i = 1;
for lambda = min_lambda : prec : max_lambda
    [x_star_i, a] = linprog(c + lambda * t, [], [], N, b, zeros(1,7), [],[],options); 
    phi_lambda(i) = a - lambda * T;
    i = i + 1;
end
plot(dom_lambda, phi_lambda);