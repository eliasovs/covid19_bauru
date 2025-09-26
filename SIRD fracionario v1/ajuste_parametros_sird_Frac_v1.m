% Elias Oliveira Vieira dos Santos
% Doutorando em Biometria, Unesp Botucatu
% elias.ov.santos@unesp.br

% O presente código 
% ajuste_parametros_sird_Frac_v1.m
% é uma ferramenta de análise epidemiológica, 
% que a partir dos dados de COVID-19 de Bauru, calibra um modelo 
% epidemiológico SIRD fracionário para que ele represente da melhor forma 
% possível os dados reais da cidade em questão.
% O código usa dados observados de casos e óbitos para encontrar os valores 
% dos parâmetros do modelo (beta, lambda, gamma) que geram as curvas de simulação 
% mais próximas da realidade observada, para diversas variações arbritárias
% de alpha. (alpha = 0.99, 0.98, 0.97, ... 1.01; 1.02, 1.03...)
% Devendo rodar o código para cada valor de alpha.

%  ========      EXECUTAR ANTES O ARQUIVO:   importardados.m      ========

%  ========      PARA RODAR O PRESENTE CÓDIGO, SÃO NECESSÁRIOS OS ARQUIVOS  
%        func_residuo_frac.m  e  fde_pi12_pc.m (Garrapa, 2018)    ========

%  = A func_residuo_frac.m, realiza os seguintes passos:
%    * A partir de um chute inicial dos parâmetros (beta, lambda, gamma),
%    e recebendo a ordem fracionária (alpha) como um valor fixo;
%    * Simula o modelo SIRD fracionário, chamando fde_pi12_pc de Garrapa (2018);
%    * Calcula os casos e óbitos acumulados a partir do resultado da simulação.x
%    * Retorna a diferença (erro) entre os dados simulados e os dados reais. 
%    onde a lsqnonlin tentará minimizar esse erro.

%  = O Arquivo SIRD___frac.m roda o modelo SIRD fracionário e gera o gráfico com as 4 variáveis =


%% === AJUSTE DOS PARÂMETROS DO MODELO SIRD FRACIONÁRIO (0 < alpha < 2) ===

% Ordem fracionária
alpha = 1.01;

% Condições iniciais
S0 = 379298;
I0 = 2;
R0 = 0;
D0 = 0;
y_0 = [S0; I0; R0; D0];

% Intervalo de tempo e passo
t0 = 0;
tf = 61;
h = 1;
t = t0:h:tf;

% Chute inicial dos parâmetros
param_ini = [0.000000480095181, 0.099847993427015, 0.011498103166975];

% Lado direito do sistema
f_fun = @(t, y, param_ini) [
    -(param_ini(1)^alpha) * y(1) * y(2);
     (param_ini(1)^alpha) * y(1) * y(2) - (param_ini(2)^alpha) * y(2) - (param_ini(3)^alpha) * y(2);
     (param_ini(2)^alpha) * y(2);
     (param_ini(3)^alpha) * y(2)
];

% Estimativa da derivada inicial (y'(0))
dy0 = [
    -param_ini(1) * S0 * I0;
     param_ini(1) * S0 * I0 - param_ini(2) * I0 - param_ini(3) * I0;
     param_ini(2) * I0;
     param_ini(3) * I0
];

%dy0 = f_fun(t0, y_inicial, param_ini);


% Conjunto de condições iniciais para alpha > 1
% (2 condições por variável: valor e derivada inicial)
y_inicial = [y_0, dy0];  % 4x2


% Recorte dos dados reais
casos_reais = casos_acumulados(1:62);
obitos_reais = obitos_acumulados(1:62);

% Limites
lb = [1e-9, 0.001, 0.001];
ub = [1e-5, 0.1, 0.05];

% Ajuste com lsqnonlin
opcoes = optimoptions('lsqnonlin', 'Display', 'iter');
[param_ajustado, resnorm] = lsqnonlin(@(param) ...
    func_residuo_frac(param, alpha, t0, tf, h, y_inicial, casos_reais, obitos_reais), ...
    param_ini, lb, ub, opcoes);

% Mostrar parâmetros ajustados
fprintf('\nParâmetros ajustados (modelo fracionário):\n');
fprintf('beta   = %.15f\n', param_ajustado(1));
fprintf('lambda = %.15f\n', param_ajustado(2));
fprintf('gamma  = %.15f\n', param_ajustado(3));

% === SIMULAÇÃO COM PARÂMETROS AJUSTADOS ===
param = param_ajustado;
[t_sim, y_sim] = fde_pi12_pc([alpha, alpha, alpha, alpha], f_fun, t0, tf, y_inicial , h, param);

S = y_sim(1,:);
I = y_sim(2,:);
R = y_sim(3,:);
D = y_sim(4,:);

casos_modelo = I + R + D;
obitos_modelo = D;

% === GRÁFICO COMPARATIVO: CASOS ===
figure;
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados reais'); hold on;
plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Modelo fracionário');
xlabel('Dias desde 30/03/2020'); ylabel('Casos acumulados'); grid on;
title('Casos acumulados de Covid-19 (modelo fracionário)');
legend('Location', 'northwest');

% === GRÁFICO COMPARATIVO: ÓBITOS ===
figure;
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Óbitos reais'); hold on;
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Modelo fracionário');
xlabel('Dias desde 30/03/2020'); ylabel('Óbitos acumulados'); grid on;
title('Óbitos acumulados de Covid-19 (modelo fracionário)');
legend('Location', 'northwest');




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

[~, y_classico] = ode45(f_classico, t, y_0);

S_class = y_classico(:,1);
I_class = y_classico(:,2);
R_class = y_classico(:,3);
D_class = y_classico(:,4);

casos_classico = I_class + R_class + D_class;
obitos_classico = D_class;


%% == GRÁFICOS COMPARANDO FRACIONÁRIO E CLÁSSICO ==
    % === CASOS ACUMULADOS ===
    figure;
    hold on; 
    plot(t, casos_classico, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
    plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Modelo fracionário');
    plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
    xlabel('Dias (desde 30/03/2020)', 'FontSize', 14);
    ylabel('Casos acumulados', 'FontSize', 14);
    title('Casos acumulados de Covid-19 em Bauru', 'FontSize', 16);
    legend('Location', 'northwest', 'FontSize', 17); 
    grid on;
    
    % === GRÁFICO: ÓBITOS ACUMULADOS ===
    figure;
    hold on;
    plot(t, obitos_classico, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
    plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Modelo fracionário');
    plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumulados de Covid-19 em Bauru'); 
    xlabel('Dias desde 30/03/2020', 'FontSize', 14);
    ylabel('Óbitos acumulados', 'FontSize', 14);
    title('Mortes acumuladas de Covid-19 em Bauru', 'FontSize', 16);
    legend('Location', 'northwest', 'FontSize', 17); grid on;
