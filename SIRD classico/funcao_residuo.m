function residuos = funcao_residuo(param, dias, casos_acumulados, obitos_acumulados, S0, I0, R0, D0)

    % Parâmetros
    beta   = param(1);
    lambda = param(2);
    gamma  = param(3);
    
    % Condições iniciais
    Y0 = [S0; I0; R0; D0];

    % Simulação
    [~, Y] = ode45(@(t, Y) sird_model(t, Y, beta, lambda, gamma), dias, Y0);

    % Compartimentos
    I = Y(:, 2);
    R = Y(:, 3);
    D = Y(:, 4);

    % Casos acumulados do modelo
    casos_modelo  = I + R + D;
    obitos_modelo = D;

    % Vetor de resíduos (diferença entre modelo e dados)
    residuos = [
        casos_modelo - casos_acumulados;
        obitos_modelo - obitos_acumulados
    ];
end
