% Elias Oliveira Vieira dos Santos
% Doutorando em Biometria, Unesp Botucatu
% elias.ov.santos@unesp.br

% O presente código 
% ajuste_parametros_sird_Frac_v2.m
% é uma ferramenta de análise epidemiológica, 
% que a partir dos dados de COVID-19 de Bauru, calibra um modelo 
% epidemiológico SIRD fracionário para que ele represente da melhor forma 
% possível os dados reais da cidade em questão.
% O código usa dados observados de casos e óbitos para encontrar os valores 
% dos parâmetros do modelo (beta, lambda, gamma) e também o parâmetro fracionário alpha,
% gerando as curvas de simulação mais próximas da realidade observada, para
% um melhor valor de alpha encontrado.

%  ========      EXECUTAR ANTES O ARQUIVO:   importardados.m      ========

%  ========      PARA RODAR O PRESENTE CÓDIGO, SÃO NECESSÁRIOS OS ARQUIVOS  
%        func_residuo_frac_alpha.m  e  fde_pi12_pc.m (Garrapa, 2018)    ========

%  = A func_residuo_frac_alpha.m, realiza os seguintes passos:
%    * Recebe um chute inicial do vetor de parâmetros, que agora possui 4
%      elementos: [beta, lambda, gamma, alpha], onde o parâmetro da ordem fracionária
%      'alpha' passa a ser uma variável a ser otimizada.
%    * Simula o modelo SIRD fracionário, chamando fde_pi12_pc de Garrapa (2018);
%    * Calcula os casos e óbitos acumulados a partir do resultado da simulação.
%    * Retorna a diferença (erro) entre os dados simulados e os dados reais. 
%    onde a lsqnonlin tentará minimizar esse erro, a partir da estimação dos 4 parâmetros.

%  = O Arquivo SIRD___frac.m roda o modelo SIRD fracionário e gera o gráfico com as 4 variáveis =


%% === AJUSTE DOS PARÂMETROS DO MODELO SIRD FRACIONÁRIO (0 < alpha < 2 ajustável) ===

%clear; clc;

% Carregar dados reais
%load('dados_sird_bauru.mat', 'casos_acumulados', 'obitos_acumulados');

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

% Estimativa inicial de parâmetros: [beta, lambda, gamma, alpha]
param_ini = [0.000000493806112, 0.0999923996574032, 0.012829998150742, 1.01];

% Limites inferiores e superiores para os parâmetros
lb = [1e-9, 0.001, 0.001, 1.001];  % alpha > 1
ub = [1e-5, 0.2, 0.05, 1.5];       % alpha máximo permitido

% Vetores reais (90 dias)
casos_reais = casos_acumulados(1:62);
obitos_reais = obitos_acumulados(1:62);

% Ajuste com lsqnonlin
opcoes = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxIterations', 400, ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10);

[param_ajustado, resnorm] = lsqnonlin(@(param) ...
    func_residuo_frac_alpha(param, t0, tf, h, y0, casos_reais, obitos_reais), ...
    param_ini, lb, ub, opcoes);

% Exibir parâmetros ajustados
fprintf('\n=== Parâmetros Ajustados (Modelo Fracionário com Alpha) ===\n');
fprintf('beta   = %.15f\n', param_ajustado(1));
fprintf('lambda = %.15f\n', param_ajustado(2));
fprintf('gamma  = %.15f\n', param_ajustado(3));
fprintf('alpha  = %.5f\n',    param_ajustado(4));

% === SIMULAÇÃO FINAL COM PARÂMETROS AJUSTADOS ===
beta   = param_ajustado(1);
lambda = param_ajustado(2);
gamma  = param_ajustado(3);
alpha  = param_ajustado(4);

% Define sistema fracionário com alpha ajustado
f_fun = @(t, y, param) [
    -(param(1)^alpha) * y(1) * y(2);
     (param(1)^alpha) * y(1) * y(2) - (param(2)^alpha) * y(2) - (param(3)^alpha) * y(2);
     (param(2)^alpha) * y(2);
     (param(3)^alpha) * y(2)
];

% Derivada inicial para alpha > 1
dy0 = [
    -beta * S0 * I0;
     beta * S0 * I0 - lambda * I0 - gamma * I0;
     lambda * I0;
     gamma * I0
];
y_inicial = [y0, dy0];

% Simulação com fde_pi12_pc
[t_sim, y_sim] = fde_pi12_pc([alpha, alpha, alpha, alpha], f_fun, t0, tf, y_inicial, h, [beta, lambda, gamma]);

S = y_sim(1,:);
I = y_sim(2,:);
R = y_sim(3,:);
D = y_sim(4,:);

casos_modelo  = I + R + D;
obitos_modelo = D;

%Dados de mortes acumuladas de Covid-19 em Bauru
%Curva estimada de mortes acumuladas de Covid-19 em Bauru
%Curva estimada de casos acumulados de Covid-19 em Bauru: modelo fracionário \alpha = 1.10835'
%Curva estimada de mortes acumuladas de Covid-19 em Bauru


% === GRÁFICO: CASOS ACUMULADOS ===
figure;
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru'); hold on;
plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo fracionário \alpha = 1.10835');
xlabel('Dias desde 30/03/2020'); ylabel('Casos acumulados'); grid on;
title('Casos acumulados de Covid-19 em Bauru');
legend('Location', 'northwest');

% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru'); hold on;
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumuladas de Covid-19 em Bauru: modelo fracionário \alpha = 1.10835');
xlabel('Dias desde 30/03/2020'); ylabel('Óbitos acumulados'); grid on;
title('Mortes acumuladas de Covid-19 em Bauru');
legend('Location', 'northwest');


%%

%% === SIMULAÇÃO DO MODELO CLÁSSICO (alpha = 1) ===
beta_cl = 0.000000480095181;
lambda_cl = 0.099847993427015;
gamma_cl = 0.011498103166975;

f_classico = @(t, y) [
    -beta_cl * y(1) * y(2);
     beta_cl * y(1) * y(2) - lambda_cl * y(2) - gamma_cl * y(2);
     lambda_cl * y(2);
     gamma_cl * y(2)
];

[~, y_classico] = ode45(f_classico, t, y0);

S_class = y_classico(:,1);
I_class = y_classico(:,2);
R_class = y_classico(:,3);
D_class = y_classico(:,4);

casos_classico = I_class + R_class + D_class;
obitos_classico = D_class;

%% === GRÁFICO: CASOS ACUMULADOS ===
figure;
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru'); hold on;
plot(t, casos_classico, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
xlabel('Dias desde 30/03/2020'); ylabel('Casos acumulados'); grid on;
title('Casos acumulados de Covid-19 em Bauru');
legend('Location', 'northwest');

%% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru'); hold on;
plot(t, obitos_classico, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumuladas de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
xlabel('Dias desde 30/03/2020'); ylabel('Óbitos acumulados'); grid on;
title('Mortes acumuladas de Covid-19 em Bauru');
legend('Location', 'northwest');

%% === GRÁFICO: CASOS ACUMULADOS ===
figure;
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru'); hold on;
plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo fracionário \alpha = 1.10835');
plot(t, casos_classico, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 14); ylabel('Casos acumulados', 'FontSize', 14); grid on;
title('Casos acumulados de Covid-19 em Bauru', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 16);

%% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru'); hold on;
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumuladas de Covid-19 em Bauru: modelo fracionário \alpha = 1.10835');
plot(t, obitos_classico, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumuladas de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 14); ylabel('Óbitos acumulados', 'FontSize', 14); grid on;
title('Mortes acumuladas de Covid-19 em Bauru', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 16);

%%
% Simulação até o dia 90, sem alterar código anterior

%tf90 = 90;
%t90 = t0:h:tf90;

%[t_sim_90, y_sim_90] = fde_pi12_pc([alpha, alpha, alpha, alpha], f_fun, t0, tf90, y_inicial, h, [beta, lambda, gamma]);

%S90 = y_sim_90(1,:);
%I90 = y_sim_90(2,:);
%R90 = y_sim_90(3,:);
%D90 = y_sim_90(4,:);

%casos_modelo_90  = I90 + R90 + D90;
%obitos_modelo_90 = D90;

% Mostrar valores no dia 90
%fprintf('\n=== Valores do modelo fracionário no dia 90 ===\n');
%fprintf('S(90) = %.2f\n', S90(end));
%fprintf('I(90) = %.2f\n', I90(end));
%fprintf('R(90) = %.2f\n', R90(end));
%fprintf('D(90) = %.2f\n', D90(end));

%fprintf('\nCasos acumulados (I+R+D) no dia 90 = %.2f\n', casos_modelo_90(end));
%fprintf('Óbitos acumulados no dia 90         = %.2f\n', obitos_modelo_90(end));



%%

%% === SIMULAÇÃO DO MODELO FRACIONÁRIO COM alpha = 1.01 ===
alpha_101 = 1.01;
beta_101 = 0.000000493806112;
lambda_101 = 0.099992399657403;
gamma_101 = 0.012829998150742;

% Sistema fracionário com alpha = 1.01
f_fun_101 = @(t, y, param) [
    -(param(1)^alpha_101) * y(1) * y(2);
     (param(1)^alpha_101) * y(1) * y(2) - (param(2)^alpha_101) * y(2) - (param(3)^alpha_101) * y(2);
     (param(2)^alpha_101) * y(2);
     (param(3)^alpha_101) * y(2)
];

% Derivadas iniciais para alpha = 1.01
dy0_101 = [
    -beta_101 * S0 * I0;
     beta_101 * S0 * I0 - lambda_101 * I0 - gamma_101 * I0;
     lambda_101 * I0;
     gamma_101 * I0
];
y_inicial_101 = [y0, dy0_101];

% Simulação
[t_101, y_101] = fde_pi12_pc([alpha_101, alpha_101, alpha_101, alpha_101], f_fun_101, t0, tf, y_inicial_101, h, [beta_101, lambda_101, gamma_101]);

% Resultados
S_101 = y_101(1,:);
I_101 = y_101(2,:);
R_101 = y_101(3,:);
D_101 = y_101(4,:);

casos_101  = I_101 + R_101 + D_101;
obitos_101 = D_101;


%%



%% === GRÁFICO: CASOS ACUMULADOS ===
figure;
hold on; % <-- aqui!
plot(t, casos_classico, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
plot(t, casos_101, 'c-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 1 \alpha = 1.01');
plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2 \alpha  = 1.10835');
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 14);
ylabel('Casos acumulados', 'FontSize', 14);
title('Casos acumulados de Covid-19 em Bauru', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 17); 
grid on;


%% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
hold on;
plot(t, obitos_classico, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
plot(t, obitos_101, 'y-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 1 \alpha = 1.01');
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2 \alpha  = 1.10835');
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumulados de Covid-19 em Bauru'); 
xlabel('Dias desde 30/03/2020', 'FontSize', 14);
ylabel('Óbitos acumulados', 'FontSize', 14);
title('Mortes acumuladas de Covid-19 em Bauru', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 17); grid on;



%% Mais gráficos


% === GRÁFICO: CASOS ACUMULADOS ===
figure;
hold on; % <-- aqui!
plot(t, casos_classico, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
%plot(t, casos_101, 'c-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 1 \alpha = 1.01');
plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2 \alpha  = 1.10835');
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11);
ylabel('Casos acumulados', 'FontSize', 11);
title('Casos acumulados de Covid-19 em Bauru', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12); 
ylim([0 400]);  % <-- Aqui você fixa os limites do eixo Y
grid on;

%%
% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
hold on;
plot(t, obitos_classico, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
%plot(t, obitos_101, 'y-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 1 \alpha = 1.01');
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2 \alpha  = 1.10835');
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumulados de Covid-19 em Bauru'); 
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11);
ylabel('Óbitos acumulados', 'FontSize', 11);
title('Mortes acumuladas de Covid-19 em Bauru', 'FontSize', 14);
ylim([0 30]);  % <-- Aqui você fixa os limites do eixo Y
legend('Location', 'northwest', 'FontSize', 12); grid on;

