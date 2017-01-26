
function [x_star, tags] = camino_mas_corto_con_cota_de_tiempo(arcos, b, c, t, T)
% CAMINO_MAS_CORTO_CON_COTA_DE_TIEMPO Esta funci?n resuelve por
% cutting-planes el problema del camino m?s corto con cotas de tiempo
% arcos = cell array con la lista de arcos del grafo, de la forma
%   inicio-fin.
% b = arreglo indicando el flujo (1 en el de salida, -1 en el de llegada, 0
%   en el resto)
% c = arreglo de costos de cada arco del grafo
% t = arreglo de tiempos de cada arco del grafo
% T = cota de tiempo
    
    % Precondiciones
    assert(length(arcos)==length(c) && length(arcos)==length(c) && length(arcos)==length(t), 'c y t deben tener la misma cantidad de elementos que arcos');
    
    num_nodos = length(b);
    
    % Si todos los elementos de mi arreglo de tiempos son mayores a T, no
    % hay soluci?n
    if all(T < t)

        % Tiro un error
        error('No existe solucion porque T es menor a cualquier t');
        
    else

        % Inicializo mi matriz nodo/arco de #nodos x #arcos
        N = zeros(num_nodos, length(arcos));

        % Por cada arco, cargo una entrada 
        for i = 1 : length(arcos)
            % me quedo con el arco actual
            arco_actual = arcos{i};
            % lo separo con el gui?n
            nodos_del_arco = strsplit(arco_actual, '-');
            % agrego las entradas correspondientes en la matriz nodo/arco
            N(str2num(nodos_del_arco{1}), i) = 1;
            N(str2num(nodos_del_arco{2}), i) = -1;
        end
        
        % Resuelvo por cutting planes
        [hubo_solucion, x_star] = resolver_por_cutting_planes(N, b, c, t, T);
        
        % si hubo solucion...
        if hubo_solucion
            tags = arcos(x_star);
        else
            tags = '';
        end
        
    end

end


function [hubo_solucion, x_star] = resolver_por_cutting_planes(N, b, c, t, T)

    % Con esta sentencia me aseguro de usar el metodo simplex, que como recorre
    % los vertices, en caso de que haya 2 caminos optimos iguales, va a
    % devolvernos uno solo. Si no usara esto, obtendria flujos fraccionarios.
    options = optimoptions('linprog', 'Algorithm', 'Simplex', 'Display', 'off');

    % Resuelvo por cutting planes
    
    lambdas = [];           %
    f_lambdas = [];         %
    phi_lambdas = [];    % estructuras para acumular valores y graficar
    xis = [];               %
    x_stars = [];           %

    % Mi epsilon para marcar convergencia
    epsilon = 1e-4;

    A = [];
    B = [];
    exit_flag = -3;
    iter = 0

    % Mientras el cutting planes nos de unbounded
    while (exit_flag==-3)

        lambda_0 = rand;
        lambda_i = lambda_0;

        % precalculo los valores de x_star, xi y f_lambda para inicializar
        lambdas = [lambdas, lambda_0];

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
        A = cat(1, A, [-1 -xi_i]); 
        % Acumulo lo que necesito para las desigualdades que resuelvo por cutting
        % planes
        B = cat(2, B, [phi_lambda_i - xi_i * lambda_i]);

        %CUTTING PLANES
        [p, f_lambda_i, exit_flag] = linprog([1 0], A, B, [], [], [-Inf 0], [], [], options);
        f_lambda_i = f_lambda_i * -1;
        lambda_i = p(end);  
        lambda_i_prev = lambda_i;

        iter = iter + 1;
    end


    % MODELO ORIGINAL
    [x_star_i, phi_lambda_i_parcial] = linprog(c + lambda_i * t, [], [], N, b, zeros(1,size(N,2)),[],[],options);
    phi_lambda_i = phi_lambda_i_parcial - lambda_i * T;

    xi_i = t * x_star_i - T;

    f_lambdas = [f_lambdas f_lambda_i];                   %
    phi_lambdas = [phi_lambdas phi_lambda_i ];   %
    xis = [xis xi_i];                                     % estruct. para 
    x_stars = [x_stars; x_star_i];                        % acumular
    lambdas = [lambdas lambda_i];                         %

    A = cat(1, A, [-1 -xi_i]);
    B = cat(2, B, phi_lambda_i - xi_i * lambda_i);

    [ phi_lambda_i,   f_lambda_i ];

    iter = iter + 1;

    while ~(abs(phi_lambda_i - f_lambda_i) < epsilon) && (iter < 30)


        %CUTTING PLANES
        [p, f_lambda_i] = linprog([1 0], A, B, [], [], [-Inf 0], [], [], options);
        f_lambda_i = f_lambda_i * -1;
        lambda_i = p(end);  
        lambda_i_prev = lambda_i;

        % MODELO ORIGINAL
        [x_star_i, phi_lambda_i_parcial] = linprog(c + lambda_i * t, [], [], N, b, zeros(1,size(N,2)),[],[],options);
        phi_lambda_i = phi_lambda_i_parcial - lambda_i * T;

        xi_i = t * x_star_i - T;

        f_lambdas = [f_lambdas f_lambda_i];                   %
        phi_lambdas = [phi_lambdas phi_lambda_i ];   %
        xis = [xis xi_i];                                     % estruct. para 
        x_stars = [x_stars; x_star_i];                        % acumular
        lambdas = [lambdas lambda_i];                         %

        A = cat(1, A, [-1 -xi_i]);
        B = cat(2, B, phi_lambda_i - xi_i * lambda_i);

        [ phi_lambda_i,   f_lambda_i ];

        iter = iter + 1;
    end

    hubo_solucion = abs(phi_lambda_i - f_lambda_i) <= epsilon;
    
    if hubo_solucion
        
        %x_star = keys(find(x_star_i));
        x_star = logical(x_star_i);

        figure
        plot(f_lambdas, 'LineWidth', 2); 
        hold on
        plot(phi_lambdas, 'LineWidth', 2)
        legend({'$f(\lambda)$', '$\phi(\lambda)$'}, 'location', 'southeast', 'Interpreter', 'LaTex');
        xlabel('iteracion', 'Interpreter', 'LaTex');
        ylabel('$f$', 'Interpreter', 'LaTex');
        xlim([1 iter]);
        grid on
        title('Evolucion del valor de $f(\lambda)$ y $\phi(\lambda)$ por iteracion', 'Interpreter', 'LaTex');
        set(findall(gcf,'-property','FontSize'),'FontSize',16)

        figure
        plot(xis, 'LineWidth', 2); 
        xlim([1 iter]);
        grid on
        title('Evolucion de $\xi$ por iteracion', 'Interpreter', 'LaTex');
        xlabel('iteracion', 'Interpreter', 'LaTex');
        ylabel('$\xi$', 'Interpreter', 'LaTex');
        set(findall(gcf,'-property','FontSize'),'FontSize',16)

        figure
        bar(x_star_i); 
        title('Solucion', 'Interpreter', 'LaTex')
        xlabel('Columna en matriz nodo/arco', 'Interpreter', 'LaTex')
        ylabel('$x^*$', 'Interpreter', 'LaTex');
        set(findall(gcf,'-property','FontSize'),'FontSize',16)

        figure
        plot(lambdas, 'LineWidth', 2);
        xlim([1 iter]);
        grid on
        title('Evolucion de $\lambda$ por iteracion', 'Interpreter', 'LaTex');
        xlabel('iteracion', 'Interpreter', 'LaTex');
        ylabel('$\lambda$', 'Interpreter', 'LaTex');
        set(findall(gcf,'-property','FontSize'),'FontSize',16)

    else
       
        x_star = [];
        
    end

    min_lambda = 0;
    prec = 0.1;
    max_lambda = 10;

    dom_lambda = min_lambda : prec : max_lambda;
    phi_lambda = zeros(length(dom_lambda), 1);
    i = 1;
    for lambda = min_lambda : prec : max_lambda
        [x_star_i, a] = linprog(c + lambda * t, [], [], N, b, zeros(1, size(N,2)), [],[],options); 
        phi_lambda(i) = a - lambda * T;
        i = i + 1;
    end
    plot(dom_lambda, phi_lambda, 'LineWidth', 2);
    xlabel('$\lambda$', 'Interpreter', 'LaTex');
    ylabel('$\phi(\lambda)$', 'Interpreter', 'LaTex');
    box on
    grid on
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
end