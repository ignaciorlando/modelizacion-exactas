
function [x_star, xi, f_lambda] = oraculo(lambda, C, t, T, N, b)
% INPUTS:
%   lambda = valor de lambda que define la funci�n objetivo
%   C = costos de cada arco
%   t = tiempo de cada arco
%   T = cotas de tiempo
%   A = matriz nodo-arco
%   b = demanda o suministro de cada nodo
% OUTPUT:
%   x_star = argmin_x ((C + \lambda * t)x - \lambda * T) 
%               s.t. N * x <= b
%   xi = subgradiente de la funci�n evaluada en x_star

    % calculo x_star
    [x_star, f_lambda] = linprog( C + lambda * t, N, b); 
    % necesito tambi�n el valor de la funci�n en x*, para eso tengo que
    % restarle el pedazo de la funci�n que no puse en el linprog
    f_lambda = f_lambda - lambda * T;
    
    % calculo xi
    xi = t * x_star + T;

end