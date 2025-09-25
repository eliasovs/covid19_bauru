function residuos = func_residuo_frac(param, alpha, t0, tf, h, y0, casos_reais, obitos_reais)

    % Função fracionária SIRD com parâmetros elevados a alpha
    f_fun = @(t,y,param) [ ...
        -(param(1)^alpha)*y(1)*y(2); ...
         (param(1)^alpha)*y(1)*y(2) - (param(2)^alpha)*y(2) - (param(3)^alpha)*y(2); ...
         (param(2)^alpha)*y(2); ...
         (param(3)^alpha)*y(2) ];

    % Resolver EDO fracionária
    [~, y] = fde_pi12_pc([alpha, alpha, alpha, alpha], f_fun, t0, tf, y0, h, param);

    % Extrair variáveis
    I = y(2, :);
    R = y(3, :);
    D = y(4, :);

    % Casos e óbitos do modelo
    casos_modelo = I + R + D;
    obitos_modelo = D;

    % Vetor de resíduos
    residuos = [
        casos_modelo(:) - casos_reais(:);
        obitos_modelo(:) - obitos_reais(:)
    ];
end
