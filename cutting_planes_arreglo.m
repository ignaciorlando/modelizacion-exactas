close all; clear; clc;
w = warning ('off','all');
epsilon = 1e-4;
% Con esta sentencia me aseguro de usar el m?todo simplex, que como recorre
% los v?rtices, en caso de que haya 2 caminos ?ptimos iguales, va a
% devolvernos uno solo. Si no usara esto, obtendr?a flujos fraccionarios.
options = optimoptions('linprog', 'Algorithm', 'Simplex', 'Display', 'off');

% PROBLEMA DE CAMINO MAS CORTO CON COTAS DE TIEMPO

% min cX
%   s.t. Nx = b
%        tx <= T
%        x >= 0

%% Defino mi grafo

% Estas son las claves de los arcos de mi grafo
keys = {'1-2', '1-3', '2-4', '2-6', '3-5', '4-6', '5-6'};

% N: Matriz nodo arco
 N = [ 1  1  0  0  0  0  0;
      -1  0  1  1  0  0  0;
       0 -1  0  0  1  0  0;
       0  0 -1  0  0  1  0;
       0  0  0  0 -1  0  1;
       0  0  0 -1  0 -1 -1; ]
   
% b: Flujo que entra y flujo que sale, indica los nodos de entrada/salida
b = [1 0 0 0 0 -1];

% c: Costos de recorrer los arcos de mi grafo
c = [2 1 2 5 2 1 2];

% t: Tiempo que lleva recorrer los arcos de mi grafo
t = [3 1 3 1 3 3 5];

% t: Tiempo m?ximo que voy a tolerar
T = 10;


%% Resoluci?n por cutting-planes

% Inicializo una estructura en la que ir guardando mis resultados parciales
to_accumulate = struct();

to_accumulate.f_lambdas = []; % llamo f_lambda a la aproximada por cutting planes
to_accumulate.x_stars = []; % son las soluciones que voy obteniendo a partir de resolver mi cutting planes
to_accumulate.xis = []; % son todos mis subgradientes

to_accumulate.lambdas = []; % para acumular los lambdas de cada iteracion
to_accumulate.phi_lambdas = []; % son mis valores reales de la funci?n objetivo en cada lambda

% mi lambda inicial es 1
lambda_i = 1;
% mi f_lambda_i es cualquiera al principio
f_lambda_i = 10;

% inicialmente A y B son vacios
A = [];
B = [];

% iter empieza en 1
iter = 1;

% Siempre entro al menos 1 vez
while (iter==1) || ~(abs(phi_lambda_i - f_lambda_i) < epsilon) 
    
    % precalculo los valores de x_star, xi y f_lambda para inicializar
    to_accumulate.lambdas = [to_accumulate.lambdas, lambda_i];

    % resuelvo y obtengo una primera solucion x* y su \phi(\lambda):
    %
    %       min    c' x + \lambda * (tx - T)
    %       x>=0
    %       N x* = b                           <- donde empieza y termina
    %
    [x_star_i, phi_lambda_i] = linprog(c + lambda_i * t, [], [], N, b, zeros(1,size(N,2)), [],[],options); 
    phi_lambda_i = phi_lambda_i - lambda_i * T;

    % Calculo el subgradiente
    xi_i = t * x_star_i - T;

    % Acumulo los subgradientes
    A = cat(1, A, [-1 -xi_i]);
    % Acumulo lo que necesito para las desigualdades que resuelvo por cutting
    % planes
    B = cat(2, B, f_lambda_i - xi_i * lambda_i);


    % ahora resuelvo por cutting planes:
    %
    %       max             phi_lambda_i
    %       \lambda >= 0
    %
    % que la unica manera que tengo de resolverlo es planteando un vector
    % p = (Z, \lambda) donde Z es una variables auxiliar y \lambda no esta en
    % mi funcion objetivo:
    %
    %       min             (1,0)' (Z, \lambda)
    %       Z, \lambda
    %       
    %       -Z - \lambda \xi_i <= \phi(lambda_i) - \psi_i \lambda_i
    %
    [p, f_lambda_i] = linprog([1 0], A, B, [], [], [-Inf 0], [], [], options);
    lambda_i = p(end);  


    % Concateno todo para ir llevando un hist?rico
    to_accumulate.f_lambdas = [to_accumulate.f_lambdas, f_lambda_i]; % llamo f_lambda a la aproximada por cutting planes
    to_accumulate.x_stars = [to_accumulate.x_stars, x_star_i]; % son las soluciones que voy obteniendo a partir de resolver mi cutting planes
    to_accumulate.xis = [to_accumulate.xis, xi_i]; % son todos mis subgradientes
    to_accumulate.phi_lambdas = [to_accumulate.phi_lambdas, phi_lambda_i]; % son mis valores reales de la funci?n objetivo en cada lambda
    
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