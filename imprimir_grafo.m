
function imprimir_grafo(num_nodos, arcos)

    % Inicializo la matriz de transicion de estados
    UG = sparse(num_nodos, num_nodos);
    
    % Por cada arco, cargo una entrada 
    for i = 1 : length(arcos)
        % me quedo con el arco actual
        arco_actual = arcos{i};
        % lo separo con el gui?n
        nodos_del_arco = strsplit(arco_actual, '-');
        % agrego las entradas correspondientes en la matriz nodo/arco
        nodo1 = str2num(nodos_del_arco{1});
        nodo2 = str2num(nodos_del_arco{2});
        UG(nodo1, nodo2) = 1;
    end

    %view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))

end