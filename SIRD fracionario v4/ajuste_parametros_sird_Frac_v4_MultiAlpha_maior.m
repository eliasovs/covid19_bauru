% Elias Oliveira Vieira dos Santos
% Doutorando em Biometria, Unesp Botucatu
% elias.ov.santos@unesp.br

% O presente código 
% ajuste_parametros_sird_Frac_v4_MultiAlpha.m
% (Alterando limites superiores dos parâmetros a serem estimados)
% é uma ferramenta de análise epidemiológica, 
% que a partir dos dados de COVID-19 de Bauru, calibra um modelo 
% epidemiológico SIRD fracionário (COM 4 ORDENS ALPHA INDEPENDENTES), para
% que ele represente da melhor forma possível os dados reais da cidade em questão.
% O código usa dados observados de casos e óbitos para encontrar os valores 
% dos parâmetros do modelo (beta, lambda, gamma) e também o parâmetros fracionários
% alpha1, alpha2, alpha3 e alpha4,
% gerando as curvas de simulação mais próximas da realidade observada.

%  ========      EXECUTAR ANTES O ARQUIVO:   importardados.m      ========

%  ========      PARA RODAR O PRESENTE CÓDIGO, SÃO NECESSÁRIOS OS ARQUIVOS  
%        func_residuo_frac_multi_alpha.m  e  fde_pi12_pc.m (Garrapa, 2018)    ========

%  = A func_residuo_frac_alpha.m, realiza os seguintes passos:
%    * Recebe um chute inicial do vetor de parâmetros, que agora possui 7
%      elementos: [beta, lambda, gamma, alpha1, alpha2, alpha3, alpha4].
%    * Nesta versão, cada equação do modelo pode ter sua própria ordem
%      fracionária, permitindo capturar diferentes 'efeitos de memória'
%      para a transmissão, recuperação e mortalidade da doença.
%    * Simula o modelo SIRD fracionário, que agora utiliza uma ordem 'alpha'
%      específica para cada equação do sistema, chamando fde_pi12_pc 
%      de Garrapa (2018);
%    * Calcula os casos e óbitos acumulados a partir do resultado da simulação.
%    * Retorna a diferença (erro) entre os dados simulados e os dados reais. 
%    onde a lsqnonlin tentará minimizar esse erro, a partir da estimação dos 7 parâmetros.

%  = O Arquivo SIRD___frac.m roda o modelo SIRD fracionário e gera o gráfico com as 4 variáveis =

%% === AJUSTE DOS PARÂMETROS DO MODELO SIRD FRACIONÁRIO COM 4 ORDENS ALPHA INDEPENDENTES ===

% Intervalo de tempo e passo
t0 = 0;
tf = 61;
h = 1;
t = t0:h:tf;

% Condições iniciais
S0 = 379298;
I0 = 2;
R0 = 0;
D0 = 0;
y0 = [S0; I0; R0; D0];

% Estimativa inicial: [beta, lambda, gamma, alpha1, alpha2, alpha3, alpha4]
%param_ini = [0.0000004938, 0.1, 0.0128, 1.1, 1.1, 1.1, 1.1];
param_ini = [0.0000004938, 0.1, 0.0128, 1.1, 1.1, 1.1, 1.1];

% Limites inferior e superior
%Estrategia 3:  
%lb = [1e-9, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
%Estrategia 4: 
lb = [1e-9, 0.001, 0.001, 1.001, 1.001, 1.001, 1.001];
ub = [1e-5, 0.2, 0.05, 1.5, 1.5, 1.5, 1.5];
%Testar com limite superior próximo de 2:
%ub = [1e-5, 0.2, 0.05, 1.99, 1.99, 1.99, 1.99];

% Dados reais
casos_reais  = casos_acumulados(1:62);
obitos_reais = obitos_acumulados(1:62);

% Opções da otimização
opcoes = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxIterations', 400, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

% Ajuste
[param_ajustado, resnorm] = lsqnonlin(@(param) func_residuo_frac_multi_alpha(param, t0, tf, h, y0, casos_reais, obitos_reais), param_ini, lb, ub, opcoes);

% Exibição
fprintf('\n=== Parâmetros Ajustados ===\n');
fprintf('beta   = %.15f\n', param_ajustado(1));
fprintf('lambda = %.15f\n', param_ajustado(2));
fprintf('gamma  = %.15f\n', param_ajustado(3));
fprintf('alpha1 = %.5f\n', param_ajustado(4));
fprintf('alpha2 = %.5f\n', param_ajustado(5));
fprintf('alpha3 = %.5f\n', param_ajustado(6));
fprintf('alpha4 = %.5f\n', param_ajustado(7));

% Parâmetros finais
beta = param_ajustado(1);
lambda = param_ajustado(2);
gamma = param_ajustado(3);
alpha1 = param_ajustado(4);
alpha2 = param_ajustado(5);
alpha3 = param_ajustado(6);
alpha4 = param_ajustado(7);

% Sistema
f_fun = @(t, y, p) [
    -(p(1)^alpha1) * y(1) * y(2);
     (p(1)^alpha2) * y(1) * y(2) - (p(2)^alpha2) * y(2) - (p(3)^alpha2) * y(2);
     (p(2)^alpha4) * y(2);
     (p(3)^alpha3) * y(2)
];

% Derivada inicial
dy0 = [
    -beta * S0 * I0;
     beta * S0 * I0 - lambda * I0 - gamma * I0;
     lambda * I0;
     gamma * I0
];
y_inicial = [y0, dy0];

% Simulação
[t_sim, y_sim] = fde_pi12_pc([alpha1, alpha2, alpha3, alpha4], f_fun, t0, tf, y_inicial, h, [beta, lambda, gamma]);
S = y_sim(1,:); I = y_sim(2,:); R = y_sim(3,:); D = y_sim(4,:);
casos_modelo = I + R + D; obitos_modelo = D;

% Gráficos
figure;
plot(t, casos_reais, 'ko', 'DisplayName', 'Dados reais'); hold on;
plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Modelo fracionário com múltiplos alphas');
xlabel('Dias'); ylabel('Casos acumulados'); title('Casos acumulados - Modelo SIRD Fracionário');
legend; grid on;

figure;
plot(t, obitos_reais, 'ko', 'DisplayName', 'Óbitos reais'); hold on;
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Modelo fracionário com múltiplos alphas');
xlabel('Dias'); ylabel('Óbitos acumulados'); title('Óbitos acumulados - Modelo SIRD Fracionário');
legend; grid on;
