function res = func_residuo_frac_alpha(param, t0, tf, h, y0, casos_reais, obitos_reais)

    % Separa os parâmetros
    beta   = param(1);
    lambda = param(2);
    gamma  = param(3);
    alpha  = param(4);

    % Define o lado direito do sistema SIRD
    f_fun = @(t, y, param_local) [
        -(param_local(1)^alpha) * y(1) * y(2);
         (param_local(1)^alpha) * y(1) * y(2) - (param_local(2)^alpha) * y(2) - (param_local(3)^alpha) * y(2);
         (param_local(2)^alpha) * y(2);
         (param_local(3)^alpha) * y(2)
    ];

    % Recalcula a derivada inicial dy0 com os novos parâmetros
    S0 = y0(1);
    I0 = y0(2);
    dy0 = [
        -beta * S0 * I0;
         beta * S0 * I0 - lambda * I0 - gamma * I0;
         lambda * I0;
         gamma * I0
    ];

    % Conjunto de condições iniciais: y(0) e y'(0)
    y_inicial = [y0, dy0];  % 4x2 matriz

    % Simulação usando fde_pi12_pc
    [~, y] = fde_pi12_pc([alpha, alpha, alpha, alpha], f_fun, t0, tf, y_inicial, h, [beta, lambda, gamma]);

    % Recupera variáveis do modelo
    I = y(2,:);
    R = y(3,:);
    D = y(4,:);

    % Casos e óbitos simulados
    casos_modelo  = I + R + D;
    obitos_modelo = D;

    % Define o vetor de resíduos: diferença entre modelo e dados
    res = [
        casos_modelo(:)  - casos_reais(:);
        obitos_modelo(:) - obitos_reais(:)
    ];
end
