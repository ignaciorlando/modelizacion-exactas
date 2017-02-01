close all;
% Este ejercicio

% Arreglo de arcos
arcos = {'9-1', '9-2', '1-3', '1-4', '2-3', '2-5', '3-4', '3-5', '4-6', '4-7', '5-6', '5-8', '6-7', '6-8', '7-10', '8-10'};
% Arreglo de costos
c = [2, 3, 1, 3, 2, 3, 3, 1, 2, 3, 2, 9, 1, 3, 12, 3];
% Arreglo de tiempos
t = [2, 2, 1, 3, 2, 3, 2, 2, 3, 3, 1, 6, 1, 5, 1, 3];
% Flujo
b = [0, 0, 0, 0, 0, 0, 0, 0, 1, -1];
% Cota de tiempo
T = 8;

imprimir_grafo(length(b), arcos);

% Resuelvo el problema
[x_star, tags] = camino_mas_corto_con_cota_de_tiempo(arcos, b, c, t, T);

% Muestro la solucion por pantalla
if isempty(x_star)
    fprintf('No hay solucion\n');
else
    fprintf('La solucion es:\n');
    disp(tags);
end