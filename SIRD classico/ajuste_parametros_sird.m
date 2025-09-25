% Elias Oliveira Vieira dos Santos
% Doutorando em Biometria, Unesp Botucatu
% elias.ov.santos@unesp.br

% O presente código é uma ferramenta de análise epidemiológica, 
% que a partir dos dados de COVID-19 de Bauru, calibra um modelo 
% epidemiológico SIRD para que ele represente da melhor forma possível os
% dados reais da cidade em questão.
% O código usa dados observados de casos e óbitos para encontrar os valores 
% dos parâmetros do modelo (beta, lambda, gamma) que geram as curvas de simulação 
% mais próximas da realidade observada.

%  ========      EXECUTAR ANTES O ARQUIVO:   importardados.m      ========

%  ========      PARA RODAR O PRESENTE CÓDIGO, SÃO NECESSÁRIOS OS ARQUIVOS  
%          sird_model.m  e funcao_residuo.m                      ========

%  = O Arquivo simuacao_sird.m roda o modelo SIRD e gera o gráfico com as 4 variáveis =

%  === SÃO NECESSÁRIOS OS ARQUIVOS  sird_model.m  e funcao_residuo.m  ===

%  = A função sird_model.m, é chamada tanto em:  
%  ajuste_parametros_sird.m, quanto na  funcao_residuo.m  para resolver
%  as equações diferenciais ordinárias clássicas (EDOs) com o ode45.   

%  = A funcao_residuo.m, realiza os seguintes passos:
%    * A partir de um chute inicial dos parâmetros (beta, lambda, gamma),
%    * Simula o modelo SIRD (chamando ode45 e sird_model);
%    * Calcula os casos e óbitos acumulados a partir do resultado da simulação.x
%    * Retorna a diferença (erro) entre os dados simulados e os dados reais. 
%    onde a lsqnonlin tentará minimizar esse erro.


%% === AJUSTE DOS PARÂMETROS DO MODELO SIRD COM LSQNONLIN ===

% Condições iniciais
S0 = 379298;
I0 = 2;
R0 = 0;
D0 = 0;

% Vetor de tempo
%dias = 0:(length(casos_acumulados)-1);
dias = 0:61;  % Usar os 62 primeiros dias

casos_acumulados  = casos_acumulados(1:62);
obitos_acumulados = obitos_acumulados(1:62);



% Chute inicial dos parâmetros [beta, lambda, gamma]
param_ini = [1e-7, 0.03, 0.005];

% Limites inferiores e superiores
lb = [1e-9, 0.001, 0.001];
ub = [1e-5, 0.1, 0.05];

% Ajuste com lsqnonlin
opcoes = optimoptions('lsqnonlin','Display','iter');
[param_ajustado, resnorm] = lsqnonlin(@(param) ...
    funcao_residuo(param, dias, casos_acumulados, obitos_acumulados, S0, I0, R0, D0), ...
    param_ini, lb, ub, opcoes);

% Mostrar parâmetros ajustados
fprintf('Parâmetros ajustados:\n');
fprintf('beta   = %.15f\n', param_ajustado(1));
fprintf('lambda = %.15f\n', param_ajustado(2));
fprintf('gamma  = %.15f\n', param_ajustado(3));


%beta   = 0.000000294858301579693;    % Taxa de contato
%lambda = 0.0345122943418598;         % Taxa de recuperação
%gamma  = 0.00835692416616626;        % Taxa de mortalidade

%param_ajustado(1) = beta
%param_ajustado(2) = lambda
%param_ajustado(3) = gamma

% Mostrar parâmetros ajustados
%fprintf('Parâmetros iniciais:\n');
%fprintf('beta   = %.15f\n', param_ajustado(1));
%fprintf('lambda = %.15f\n', param_ajustado(2));
%fprintf('gamma  = %.15f\n', param_ajustado(3));

%% === SIMULAÇÃO COM OS PARÂMETROS AJUSTADOS ===

% Simular com parâmetros ajustados
[~, Y_ajustado] = ode45(@(t, Y) sird_model(t, Y, ...
    param_ajustado(1), param_ajustado(2), param_ajustado(3)), dias, [S0; I0; R0; D0]);

% Extrair compartimentos
I_aj = Y_ajustado(:, 2);
R_aj = Y_ajustado(:, 3);
D_aj = Y_ajustado(:, 4);

% Casos acumulados simulados
casos_modelo_ajustado  = I_aj + R_aj + D_aj;
obitos_modelo_ajustado = D_aj;

% === GRÁFICO 1: CASOS ACUMULADOS ===
figure;
plot(dias, casos_acumulados, 'ko', 'MarkerSize', 6, ...
    'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
hold on;
plot(dias, casos_modelo_ajustado, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru');
xlabel('Dias desde 30/03/2020');
ylabel('Casos Acumulados');
title('Curva estimada de casos acumulados de Covid-19 em Bauru');
legend('Location', 'northwest');
grid on;

% === GRÁFICO 2: ÓBITOS ACUMULADOS ===
figure;
plot(dias, obitos_acumulados, 'ko', 'MarkerSize', 6, ...
    'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru');
hold on;
plot(dias, obitos_modelo_ajustado, 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Curva estimada de mortes acumuladas de Covid-19 em Bauru');
xlabel('Dias desde 30/03/2020');
ylabel('Óbitos Acumulados');
title('Curva estimada de mortes acumuladas de Covid-19 em Bauru');
legend('Location', 'northwest');
grid on;


%% === SIMULAÇÃO COM OS PARÂMETROS ANTERIORES (FORNECIDOS MANUALMENTE) ===
beta_ant   = 0.000000294858301579693;
lambda_ant = 0.0345122943418598;
gamma_ant  = 0.00835692416616626;

[~, Y_anterior] = ode45(@(t, Y) sird_model(t, Y, ...
    beta_ant, lambda_ant, gamma_ant), dias, [S0; I0; R0; D0]);

I_ant = Y_anterior(:, 2);
R_ant = Y_anterior(:, 3);
D_ant = Y_anterior(:, 4);

casos_modelo_anterior  = I_ant + R_ant + D_ant;
obitos_modelo_anterior = D_ant;

%% === GRÁFICO 1: CASOS ACUMULADOS (3 curvas) ===
figure;
hold on;
plot(dias, casos_modelo_ajustado, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru');
plot(dias, casos_modelo_anterior, 'g-', 'LineWidth', 2, ...
    'DisplayName', 'Curva de casos acumulados com parâmetros estimados em Santos (2024)');
plot(dias, casos_acumulados, 'ko', 'MarkerSize', 6, ...
    'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 14);
ylabel('Casos Acumulados', 'FontSize', 14);
title('Curva estimada de casos acumulados de Covid-19 em Bauru', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 17);
grid on;

%% === GRÁFICO 2: ÓBITOS ACUMULADOS (3 curvas) ===
figure;
hold on;
plot(dias, obitos_modelo_ajustado, 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Curva estimada de mortes acumuladas de Covid-19 em Bauru');
plot(dias, obitos_modelo_anterior, 'm-', 'LineWidth', 2, ...
    'DisplayName', 'Curva de mortes acumuladas com parâmetros estimados em Santos (2024)');
plot(dias, obitos_acumulados, 'ko', 'MarkerSize', 6, ...
    'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 14);
ylabel('Óbitos Acumulados', 'FontSize', 14);
title('Curva estimada de mortes acumuladas de Covid-19 em Bauru', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 17);
grid on;
